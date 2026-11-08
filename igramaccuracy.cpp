// Accuracy of the average of several interferograms taken in ONE direction.
//
// Every shot i gives a wave front W_i (piston, tilt and focus removed: they
// change from shot to shot).  The fringes are turned between shots, and the
// direction of the tilt phi_i that made the fringes is still in the data, so
// the fringe orientation of each shot is known without asking the user.
//
// Model, per pixel:
//     W_i = M + cos(phi_i) P1 + sin(phi_i) Q1 + cos(2 phi_i) P2 + sin(2 phi_i) Q2 + noise_i
// M is the mirror, the P/Q maps are errors that turn with the fringes.
// On real data (AR2627/55-0, AR2623/36-0, AR2612/42-0, Nov 2026) every session
// had 28-40 nm of astigmatism turning with the tilt direction (period 360
// degrees).  A set of shots spread evenly around the circle cancels it in the
// plain average; an unbalanced set does not, and the shift is reported.
//
// The noise is what is left after the model.  Consecutive shots that see the
// same air are correlated; the correlation of consecutive residuals gives an
// effective number of independent shots, N(1-r)/(1+r).  The correlation has
// to be taken after the model: the two shots of a pair share their angle and
// look correlated through the fringe term alone.
//
// Errors that are the same in every shot (reference, null, the stand, steady
// air currents) are invisible here - only measuring another way shows them.

#include "igramaccuracy.h"
#include "wavefront.h"
#include "mirrordlg.h"

#include <opencv2/opencv.hpp>
#include <QApplication>
#include <QDesktopServices>
#include <QDialog>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QObject>
#include <QPainter>
#include <QPageLayout>
#include <QPageSize>
#include <QPrinter>
#include <QPushButton>
#include <QTextBrowser>
#include <QTextDocument>
#include <QUrl>
#include <QVBoxLayout>
#include <QVariant>
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

qreal reportImageDpi();   // surfacemanager.cpp

namespace {

const int NZ = 11;
const char *zNames[NZ] = {"Piston", "Tilt X", "Tilt Y", "Focus", "Astig X", "Astig Y",
                          "Coma X", "Coma Y", "Spherical", "Trefoil X", "Trefoil Y"};

// Zernike terms 0..10 in DFTFringe's order and normalisation.
void basis(double rho, double th, double *b){
    double r2 = rho * rho;
    b[0] = 1.;
    b[1] = rho * cos(th);
    b[2] = rho * sin(th);
    b[3] = 2. * r2 - 1.;
    b[4] = r2 * cos(2. * th);
    b[5] = r2 * sin(2. * th);
    b[6] = (3. * r2 - 2.) * rho * cos(th);
    b[7] = (3. * r2 - 2.) * rho * sin(th);
    b[8] = 6. * r2 * r2 - 6. * r2 + 1.;
    b[9] = r2 * rho * cos(3. * th);
    b[10] = r2 * rho * sin(3. * th);
}

double median(std::vector<double> v){
    if (v.empty())
        return 0.;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double rowRms(const cv::Mat &m, int row){
    cv::Mat r = m.row(row);
    return sqrt(r.dot(r) / r.cols);
}

} // namespace

AccuracyResult computeIgramAccuracy(const QList<wavefront *> &wfs){
    AccuracyResult res;
    const int N = wfs.size();
    res.n = N;
    if (N < 4){
        res.error = QObject::tr("Select at least four wave fronts.");
        return res;
    }
    if (mirrorDlg::get_Instance()->isEllipse()){
        res.error = QObject::tr("The accuracy estimate works only for round mirrors.");
        return res;
    }
    const cv::Size size = wfs[0]->data.size();
    for (wavefront *wf : wfs){
        if (wf->data.size() != size){
            res.error = QObject::tr("The selected wave fronts have different sizes.");
            return res;
        }
    }

    // Pixels valid in every shot.
    cv::Mat valid(size, CV_8U, cv::Scalar(255));
    for (wavefront *wf : wfs){
        const cv::Mat &m = (wf->workMask.size() == size) ? wf->workMask : wf->mask;
        if (m.size() != size){
            res.error = QObject::tr("A selected wave front has no mask.");
            return res;
        }
        cv::bitwise_and(valid, m, valid);
    }
    const double cx = wfs[0]->m_outside.m_center.x();
    const double cy = wfs[0]->m_outside.m_center.y();
    const double rad = wfs[0]->m_outside.m_radius;
    if (rad <= 0.){
        res.error = QObject::tr("The selected wave fronts have no outline.");
        return res;
    }
    const int step = std::max(1, size.width / 200);
    std::vector<cv::Point> pts;
    std::vector<double> rows;
    for (int y = 0; y < size.height; y += step){
        for (int x = 0; x < size.width; x += step){
            if (valid.at<uchar>(y, x) == 0)
                continue;
            double ux = (x - cx) / rad, uy = (y - cy) / rad;
            double rho = sqrt(ux * ux + uy * uy);
            if (rho > 1.)
                continue;
            double b[NZ];
            basis(rho, atan2(uy, ux), b);
            pts.push_back(cv::Point(x, y));
            rows.insert(rows.end(), b, b + NZ);
        }
    }
    const int P = pts.size();
    if (P < 500){
        res.error = QObject::tr("Too few pixels are valid in all the selected wave fronts.");
        return res;
    }
    cv::Mat B(P, NZ, CV_64F, rows.data());
    cv::Mat Bpinv;
    cv::invert(B, Bpinv, cv::DECOMP_SVD);          // NZ x P

    // Each shot: Zernike terms and the map with piston, tilt and focus removed (nm).
    cv::Mat W(N, P, CV_64F);
    cv::Mat Z(N, NZ, CV_64F);
    res.lambdaNm = wfs[0]->lambda > 0. ? wfs[0]->lambda : mirrorDlg::get_Instance()->lambda;
    for (int i = 0; i < N; ++i){
        const double lam = wfs[i]->lambda > 0. ? wfs[i]->lambda : res.lambdaNm;
        cv::Mat v(P, 1, CV_64F);
        for (int k = 0; k < P; ++k)
            v.at<double>(k) = wfs[i]->data(pts[k].y, pts[k].x) * lam;
        cv::Mat z = Bpinv * v;
        cv::Mat w = v - B.colRange(0, 4) * z.rowRange(0, 4);
        cv::Mat(w.t()).copyTo(W.row(i));
        cv::Mat(z.t()).copyTo(Z.row(i));
    }

    // All shots must have the spherical term of the same sign (Average checks the same).
    int neg = 0;
    for (int i = 0; i < N; ++i)
        if (Z.at<double>(i, 8) < 0.)
            ++neg;
    if (neg != 0 && neg != N){
        res.error = QObject::tr("Some of the selected wave fronts are inverted relative to the others. "
                                "Invert them first.");
        return res;
    }

    // Tilt direction = fringe orientation.
    std::vector<double> phi(N), tilt(N);
    for (int i = 0; i < N; ++i){
        phi[i] = atan2(Z.at<double>(i, 2), Z.at<double>(i, 1));
        tilt[i] = std::hypot(Z.at<double>(i, 1), Z.at<double>(i, 2)) / res.lambdaNm;
    }
    {
        std::vector<double> a(N);
        for (int i = 0; i < N; ++i)
            a[i] = fmod(phi[i] * 180. / M_PI + 360., 360.);
        std::sort(a.begin(), a.end());
        double gap = a[0] + 360. - a[N - 1];
        for (int i = 1; i < N; ++i)
            gap = std::max(gap, a[i] - a[i - 1]);
        res.maxGapDeg = gap;
        double c1 = 0, s1 = 0, c2 = 0, s2 = 0;
        for (int i = 0; i < N; ++i){
            c1 += cos(phi[i]); s1 += sin(phi[i]);
            c2 += cos(2 * phi[i]); s2 += sin(2 * phi[i]);
        }
        res.balance1 = std::hypot(c1, s1) / N;
        res.balance2 = std::hypot(c2, s2) / N;
    }

    // The fringe terms need enough shots and angles all around.
    const double medTilt = median(tilt);
    res.fringeModel = true;
    if (N < 8){
        res.fringeModel = false;
        res.warnings << QObject::tr("Fewer than 8 shots: the errors that turn with the fringes "
                                    "cannot be separated from the noise.");
    }
    if (medTilt < 2.){
        res.fringeModel = false;
        res.warnings << QObject::tr("The wave fronts have almost no tilt left (%1 waves), so the fringe "
                                    "direction of each shot is unknown.").arg(medTilt, 0, 'f', 1);
    }
    if (res.maxGapDeg > 120.){
        res.fringeModel = false;
        res.warnings << QObject::tr("The fringe directions leave a gap of %1 degrees. Errors that turn "
                                    "with the fringes do not cancel and cannot be estimated - the result "
                                    "may be off by more than the noise figures say. Spread the fringe "
                                    "directions around the whole circle.").arg(res.maxGapDeg, 0, 'f', 0);
    }

    const int K = res.fringeModel ? 5 : 1;
    cv::Mat X(N, K, CV_64F);
    for (int i = 0; i < N; ++i){
        X.at<double>(i, 0) = 1.;
        if (K == 5){
            X.at<double>(i, 1) = cos(phi[i]);
            X.at<double>(i, 2) = sin(phi[i]);
            X.at<double>(i, 3) = cos(2 * phi[i]);
            X.at<double>(i, 4) = sin(2 * phi[i]);
        }
    }
    cv::Mat Xpinv, XtXinv;
    cv::invert(X, Xpinv, cv::DECOMP_SVD);            // K x N
    cv::invert(X.t() * X, XtXinv, cv::DECOMP_SVD);
    const int dof = N - K;

    // Plain scatter about the average.
    cv::Mat avg;
    cv::reduce(W, avg, 0, cv::REDUCE_AVG);
    double ss0 = 0.;
    for (int i = 0; i < N; ++i){
        cv::Mat d = W.row(i) - avg;
        ss0 += d.dot(d) / P;
    }
    res.sigmaPlainNm = sqrt(ss0 / (N - 1));

    // Model fit and its residuals.
    cv::Mat coef = Xpinv * W;                         // K x P
    cv::Mat R = W - X * coef;
    std::vector<double> rms(N);
    double ss = 0.;
    for (int i = 0; i < N; ++i){
        rms[i] = rowRms(R, i);
        ss += rms[i] * rms[i];
    }
    res.sigmaNm = sqrt(ss / dof);
    if (K > 1){
        cv::Mat L = X.colRange(1, K) * coef.rowRange(1, K);   // N x P, the turning part
        res.lockedNm = sqrt(L.dot(L) / (double(N) * P));
        cv::Mat b;
        cv::reduce(L, b, 0, cv::REDUCE_AVG);                  // plain average minus M
        cv::Scalar m = cv::mean(b);
        cv::Mat bc = b - m.val[0];
        res.biasNm = sqrt(bc.dot(bc) / P);
    }

    // Correlation of consecutive residuals -> effective number of shots.
    // Neighbours with nearly the same fringe direction show the air best; the
    // fit itself makes neighbours of different direction slightly negatively
    // correlated (42-0: 0.26 within pairs, 0.06 over all neighbours).
    double csSame = 0., csAll = 0.;
    int nSame = 0;
    for (int i = 0; i + 1 < N; ++i){
        double d = sqrt(R.row(i).dot(R.row(i)) * R.row(i + 1).dot(R.row(i + 1)));
        double c = d > 0. ? R.row(i).dot(R.row(i + 1)) / d : 0.;
        csAll += c;
        double dphi = fabs(std::remainder(phi[i + 1] - phi[i], 2. * M_PI));
        if (dphi < 20. * M_PI / 180.){
            csSame += c;
            ++nSame;
        }
    }
    res.corr = nSame >= 3 ? csSame / nSame : csAll / (N - 1);
    res.nEff = res.corr > 0. ? N * (1. - res.corr) / (1. + res.corr) : N;
    res.nEff = std::max(2., std::min<double>(N, res.nEff));
    res.floorNm = res.sigmaNm / sqrt(res.nEff);
    const double inflate = sqrt(N / res.nEff);

    // RMS of the average without spherical (the null is not applied here).
    {
        cv::Mat za = Bpinv * avg.t();
        cv::Mat low = B.col(0) * za.at<double>(0) + B.col(1) * za.at<double>(1) + B.col(2) * za.at<double>(2)
                    + B.col(3) * za.at<double>(3) + B.col(8) * za.at<double>(8);
        cv::Mat d = avg.t() - low;
        res.rmsAvgNm = sqrt(d.dot(d) / P);
    }

    // Shots.
    const double med = median(rms);
    for (int i = 0; i < N; ++i){
        AccuracyShot s;
        s.name = QFileInfo(wfs[i]->name).fileName();
        s.tiltDeg = fmod(phi[i] * 180. / M_PI + 360., 360.);
        s.tiltWaves = tilt[i];
        s.residualNm = rms[i];
        s.ratio = med > 0. ? rms[i] / med : 0.;
        s.flag = s.ratio >= 2.5 ? 2 : (s.ratio > 1.5 ? 1 : 0);
        for (int k = 0; k < NZ; ++k)
            s.z[k] = Z.at<double>(i, k);
        res.shots.push_back(s);
        if (s.flag == 2)
            res.warnings << QObject::tr("%1 is far out (%2 x the median) - measure it again.")
                            .arg(s.name).arg(s.ratio, 0, 'f', 1);
    }

    // Zernike terms: plain mean (what Average gives) and its standard error.
    std::vector<double> se(NZ, 0.);
    std::vector<cv::Mat> lockedTerm(NZ);
    for (int k = 4; k < NZ; ++k){
        cv::Mat y = Z.col(k).clone();
        cv::Mat c = Xpinv * y;
        cv::Mat r = y - X * c;
        double s2 = r.dot(r) / dof;
        AccuracyTerm t;
        t.name = zNames[k];
        t.mean = cv::mean(y).val[0];
        t.sigma1 = sqrt(s2);
        t.se = sqrt(s2 * XtXinv.at<double>(0, 0)) * inflate;
        se[k] = t.se;
        res.terms.push_back(t);
        if (K > 1)
            lockedTerm[k] = X.colRange(1, K) * c.rowRange(1, K);
    }
    double ax = res.terms[0].mean, ay = res.terms[1].mean;
    res.astigMean = std::hypot(ax, ay);
    res.astigCorrected = sqrt(std::max(0., ax * ax + ay * ay - se[4] * se[4] - se[5] * se[5]));
    res.astigSe = sqrt(0.5 * (se[4] * se[4] + se[5] * se[5]));
    res.astigAngleDeg = 0.5 * atan2(ay, ax) * 180. / M_PI;
    if (K > 1){
        double s = 0.;
        for (int i = 0; i < N; ++i)
            s += pow(lockedTerm[4].at<double>(i), 2) + pow(lockedTerm[5].at<double>(i), 2);
        res.astigLockedNm = sqrt(s / N);
    }

    // Leave one out.
    res.looAstigMin = res.looSphMin = 1e300;
    res.looAstigMax = res.looSphMax = -1e300;
    double sx = 0, sy = 0, ssph = 0;
    for (int i = 0; i < N; ++i){
        sx += Z.at<double>(i, 4); sy += Z.at<double>(i, 5); ssph += Z.at<double>(i, 8);
    }
    for (int i = 0; i < N; ++i){
        double a = std::hypot(sx - Z.at<double>(i, 4), sy - Z.at<double>(i, 5)) / (N - 1);
        double s = (ssph - Z.at<double>(i, 8)) / (N - 1);
        res.looAstigMin = std::min(res.looAstigMin, a);
        res.looAstigMax = std::max(res.looAstigMax, a);
        res.looSphMin = std::min(res.looSphMin, s);
        res.looSphMax = std::max(res.looSphMax, s);
    }

    if (res.nEff < 0.7 * N)
        res.warnings << QObject::tr("Consecutive shots are alike (correlation %1): slow air movement. "
                                    "They count as %2 independent shots, not %3.")
                        .arg(res.corr, 0, 'f', 2).arg(res.nEff, 0, 'f', 1).arg(N);
    if (res.fringeModel && res.biasNm > res.floorNm)
        res.warnings << QObject::tr("The fringe directions are unbalanced: the plain average is shifted "
                                    "by about %1 nm, more than the noise floor.").arg(res.biasNm, 0, 'f', 1);

    res.valid = true;
    return res;
}

namespace {

// Astigmatism of each shot, X against Y.
QImage astigPlot(const AccuracyResult &r){
    const int S = 900;
    QImage img(S, S, QImage::Format_RGB32);
    img.fill(Qt::white);
    QPainter p(&img);
    p.setRenderHint(QPainter::Antialiasing);
    double lim = 1.;
    for (const AccuracyShot &s : r.shots)
        lim = std::max(lim, std::max(fabs(s.z[4]), fabs(s.z[5])));
    lim = std::max(lim, r.astigMean + 2 * r.astigSe);
    lim *= 1.15;
    const double m = 90.;
    const double sc = (S / 2. - m) / lim;
    auto toPx = [&](double x, double y){ return QPointF(S / 2. + x * sc, S / 2. - y * sc); };
    QFont f = p.font();
    f.setPixelSize(24);
    p.setFont(f);
    p.setPen(QPen(QColor(200, 200, 200), 2));
    p.drawLine(QPointF(m, S / 2.), QPointF(S - m, S / 2.));
    p.drawLine(QPointF(S / 2., m), QPointF(S / 2., S - m));
    p.setPen(Qt::black);
    p.drawText(QRectF(S - m - 200, S / 2. + 5, 200, 30), Qt::AlignRight, QString("Astig X"));
    p.drawText(QRectF(S / 2. + 8, m - 35, 200, 30), Qt::AlignLeft, QString("Astig Y"));
    p.drawText(QRectF(m, S - m + 20, S - 2 * m, 30), Qt::AlignCenter,
               QString("each mark one shot, nm wavefront, axes to %1 nm").arg(lim, 0, 'f', 0));
    // 2 standard errors around the mean
    QPointF c = toPx(r.terms[0].mean, r.terms[1].mean);
    p.setPen(QPen(QColor(0, 0, 0), 2, Qt::DashLine));
    p.drawEllipse(c, 2 * r.astigSe * sc, 2 * r.astigSe * sc);
    for (size_t i = 0; i < r.shots.size(); ++i){
        const AccuracyShot &s = r.shots[i];
        QColor col = s.flag == 2 ? QColor(210, 30, 30) : (s.flag == 1 ? QColor(230, 140, 0) : QColor(30, 90, 200));
        QPointF q = toPx(s.z[4], s.z[5]);
        p.setPen(Qt::NoPen);
        p.setBrush(col);
        p.drawEllipse(q, 8, 8);
        p.setPen(col);
        p.drawText(q + QPointF(10, -8), QString::number(i + 1));
    }
    p.setPen(QPen(Qt::black, 4));
    p.drawLine(c + QPointF(-14, 0), c + QPointF(14, 0));
    p.drawLine(c + QPointF(0, -14), c + QPointF(0, 14));
    return img;
}

QString nm(double v, int d = 1){ return QString::number(v, 'f', d); }

} // namespace

QString accuracyReportHtml(const AccuracyResult &r, const QString &title){
    QString h;
    h += "<html><body style='font-family:Arial; font-size:10pt'>";
    h += QString("<h2>Accuracy of the average</h2><p>%1</p>").arg(title.toHtmlEscaped());
    h += QString("<p>%1 shots, wavelength %2 nm. All values are nm of <b>wavefront</b> "
                 "(surface is half), without smoothing.</p>").arg(r.n).arg(r.lambdaNm, 0, 'f', 1);

    if (!r.warnings.isEmpty()){
        h += "<p style='color:#b00000'><b>Attention</b></p><ul>";
        for (const QString &w : r.warnings)
            h += "<li style='color:#b00000'>" + w.toHtmlEscaped() + "</li>";
        h += "</ul>";
    }

    h += "<h3>Noise</h3><table border='1' cellspacing='0' cellpadding='4'>";
    h += QString("<tr><td><b>Noise floor of the average</b></td><td><b>%1 nm</b></td></tr>").arg(nm(r.floorNm));
    h += QString("<tr><td>Scatter of a single shot (random part)</td><td>%1 nm</td></tr>").arg(nm(r.sigmaNm));
    h += QString("<tr><td>Scatter of a single shot about the plain average</td><td>%1 nm</td></tr>").arg(nm(r.sigmaPlainNm));
    h += QString("<tr><td>Independent shots (consecutive correlation %1)</td><td>%2 of %3</td></tr>")
         .arg(nm(r.corr, 2)).arg(nm(r.nEff)).arg(r.n);
    if (r.fringeModel){
        h += QString("<tr><td>Error turning with the fringes, per shot</td><td>%1 nm</td></tr>").arg(nm(r.lockedNm));
        h += QString("<tr><td>Of it left in the average (unbalanced directions)</td><td>%1 nm</td></tr>").arg(nm(r.biasNm));
    }
    h += QString("<tr><td>Balance of the fringe directions |&Sigma;e<sup>i&phi;</sup>|/N, |&Sigma;e<sup>2i&phi;</sup>|/N</td>"
                 "<td>%1, %2 (0 = perfect)</td></tr>").arg(nm(r.balance1, 2)).arg(nm(r.balance2, 2));
    h += QString("<tr><td>Largest gap between fringe directions</td><td>%1&deg;</td></tr>").arg(nm(r.maxGapDeg, 0));
    h += QString("<tr><td>RMS of the average (piston, tilt, focus, spherical removed, no null)</td><td>%1 nm</td></tr>")
         .arg(nm(r.rmsAvgNm));
    h += "</table>";

    h += "<p>The noise floor is how much the average itself may still be off because of random errors "
         "(air, vibration, analysis). Errors that are the same in every shot - the stand, the reference, "
         "steady air currents - do not show here.</p>";

    h += "<h3>Zernike terms of the average</h3><table border='1' cellspacing='0' cellpadding='4'>";
    h += "<tr><th>Term</th><th>Average</th><th>&plusmn; std. error</th><th>single shot scatter</th></tr>";
    for (const AccuracyTerm &t : r.terms)
        h += QString("<tr><td>%1</td><td align='right'>%2</td><td align='right'>%3</td><td align='right'>%4</td></tr>")
             .arg(t.name).arg(nm(t.mean)).arg(nm(t.se)).arg(nm(t.sigma1));
    h += "</table>";
    h += "<p>Spherical is before the null. The std. error is 1 sigma; about 95% of repeated averages fall within "
         "2 std. errors.</p>";
    h += QString("<p><b>Astigmatism</b> %1 &plusmn; %2 nm (corrected for noise: %3 nm). Leaving out one shot it "
                 "ranges %4 &ndash; %5 nm.").arg(nm(r.astigMean)).arg(nm(r.astigSe)).arg(nm(r.astigCorrected))
         .arg(nm(r.looAstigMin)).arg(nm(r.looAstigMax));
    if (r.fringeModel)
        h += QString(" Astigmatism turning with the fringes: %1 nm per shot - it cancels when the fringe directions "
                     "are spread evenly.").arg(nm(r.astigLockedNm));
    h += QString(" Spherical leaving out one shot: %1 &ndash; %2 nm.</p>").arg(nm(r.looSphMin)).arg(nm(r.looSphMax));
    h += "<p><img src='accuracy://astig'></p>";

    h += "<h3>Shots</h3><table border='1' cellspacing='0' cellpadding='3'>";
    h += "<tr><th>#</th><th>Wave front</th><th>Fringe dir.</th><th>Tilt, waves</th><th>Deviation, nm</th>"
         "<th>&times; median</th><th></th><th>Astig X</th><th>Astig Y</th><th>Spherical</th></tr>";
    for (size_t i = 0; i < r.shots.size(); ++i){
        const AccuracyShot &s = r.shots[i];
        QString flag = s.flag == 2 ? "<b style='color:#c00000'>far out</b>"
                     : (s.flag == 1 ? "<span style='color:#c07000'>deviating</span>" : "");
        h += QString("<tr><td>%1</td><td>%2</td><td align='right'>%3&deg;</td><td align='right'>%4</td>"
                     "<td align='right'>%5</td><td align='right'>%6</td><td>%7</td><td align='right'>%8</td>"
                     "<td align='right'>%9</td><td align='right'>%10</td></tr>")
             .arg(i + 1).arg(s.name.toHtmlEscaped()).arg(nm(s.tiltDeg, 0)).arg(nm(s.tiltWaves))
             .arg(nm(s.residualNm)).arg(nm(s.ratio, 2)).arg(flag)
             .arg(nm(s.z[4])).arg(nm(s.z[5])).arg(nm(s.z[8]));
    }
    h += "</table>";
    h += "<p>Deviation is what is left of each shot after the average and the error that turns with the fringes. "
         "Deviating: over 1.5&times; the median. Far out: 2.5&times; or more - worth measuring again. "
         "No shot is left out of the average automatically.</p>";
    h += "</body></html>";
    return h;
}

QString showAccuracyReport(const AccuracyResult &r, const QString &title,
                           const QString &pdfPath, QWidget *parent){
    const QString html = accuracyReportHtml(r, title);
    const QImage plot = astigPlot(r);

    // PDF.
    QString written;
    QString pdfError;
    if (!pdfPath.isEmpty()){
        QFile test(pdfPath);
        if (test.exists() && !test.open(QIODevice::ReadWrite)){
            pdfError = QObject::tr("Could not write %1 - is it open in a PDF viewer?").arg(QDir::toNativeSeparators(pdfPath));
        }
        else {
            test.close();
            QPrinter printer(QPrinter::HighResolution);
            printer.setColorMode(QPrinter::Color);
            printer.setOutputFormat(QPrinter::PdfFormat);
            printer.setPageSize(QPageSize(QPageSize::A4));
            printer.setOutputFileName(pdfPath);
            // Image sized in the DPI QTextDocument prints with, 60% of the text width.
            double textWidthIn = printer.pageLayout().pageSize().size(QPageSize::Inch).width() - 2 * 2 / 2.54;
            QTextDocument doc;
            doc.addResource(QTextDocument::ImageResource, QUrl("accuracy://astig"),
                            QVariant(plot.scaledToWidth(qRound(0.6 * textWidthIn * reportImageDpi()),
                                                        Qt::SmoothTransformation)));
            doc.setHtml(html);
            doc.print(&printer);
            written = pdfPath;
        }
    }

    // Window.
    QDialog dlg(parent);
    dlg.setWindowTitle(QObject::tr("Accuracy of the average"));
    dlg.resize(950, 850);
    QVBoxLayout *lay = new QVBoxLayout(&dlg);
    QTextBrowser *tb = new QTextBrowser(&dlg);
    tb->document()->addResource(QTextDocument::ImageResource, QUrl("accuracy://astig"),
                                QVariant(plot.scaledToWidth(500, Qt::SmoothTransformation)));
    tb->setHtml(html);
    lay->addWidget(tb);
    QHBoxLayout *row = new QHBoxLayout;
    QLabel *where = new QLabel(&dlg);
    where->setTextInteractionFlags(Qt::TextSelectableByMouse);
    if (!written.isEmpty())
        where->setText(QObject::tr("Saved: %1").arg(QDir::toNativeSeparators(written)));
    else if (!pdfError.isEmpty())
        where->setText("<span style='color:#b00000'>" + pdfError.toHtmlEscaped() + "</span>");
    else
        where->setText(QObject::tr("No folder known for the PDF."));
    row->addWidget(where, 1);
    if (!written.isEmpty()){
        QPushButton *open = new QPushButton(QObject::tr("Open PDF"), &dlg);
        QObject::connect(open, &QPushButton::clicked, [written](){
            QDesktopServices::openUrl(QUrl::fromLocalFile(written));
        });
        row->addWidget(open);
    }
    QPushButton *close = new QPushButton(QObject::tr("Close"), &dlg);
    QObject::connect(close, &QPushButton::clicked, &dlg, &QDialog::accept);
    row->addWidget(close);
    lay->addLayout(row);
    dlg.exec();
    return written;
}
