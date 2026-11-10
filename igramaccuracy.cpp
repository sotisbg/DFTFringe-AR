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
            res.warnings << QObject::tr("%1 is far out (%2 x the median) - delete it from the wave front list "
                                        "and run Average and Accuracy again without it.")
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
    double looMove = -1.;
    for (int i = 0; i < N; ++i){
        double a = std::hypot(sx - Z.at<double>(i, 4), sy - Z.at<double>(i, 5)) / (N - 1);
        double s = (ssph - Z.at<double>(i, 8)) / (N - 1);
        res.looAstigMin = std::min(res.looAstigMin, a);
        res.looAstigMax = std::max(res.looAstigMax, a);
        res.looSphMin = std::min(res.looSphMin, s);
        res.looSphMax = std::max(res.looSphMax, s);
        if (fabs(a - res.astigMean) > looMove){
            looMove = fabs(a - res.astigMean);
            res.looAstigShot = i;
        }
    }

    // Without the far-out shots: what Average would give if they are left out.
    {
        std::vector<int> keep;
        for (int i = 0; i < N; ++i)
            if (res.shots[i].flag != 2)
                keep.push_back(i);
        const int M = keep.size();
        if (M < N && M >= K + 3){
            cv::Mat Xc(M, K, CV_64F), Zc(M, NZ, CV_64F);
            for (int j = 0; j < M; ++j){
                X.row(keep[j]).copyTo(Xc.row(j));
                Z.row(keep[j]).copyTo(Zc.row(j));
            }
            cv::Mat Xcp, XcI;
            cv::invert(Xc, Xcp, cv::DECOMP_SVD);
            cv::invert(Xc.t() * Xc, XcI, cv::DECOMP_SVD);
            for (int k = 4; k < NZ; ++k){
                cv::Mat y = Zc.col(k).clone();
                cv::Mat c = Xcp * y;
                cv::Mat r = y - Xc * c;
                double s2 = r.dot(r) / (M - K);
                res.termCleanMean.push_back(cv::mean(y).val[0]);
                res.termCleanSe.push_back(sqrt(s2 * XcI.at<double>(0, 0)) * inflate);
            }
            res.hasClean = true;
            res.nClean = M;
            res.astigClean = std::hypot(res.termCleanMean[0], res.termCleanMean[1]);
            res.astigCleanSe = sqrt(0.5 * (res.termCleanSe[0] * res.termCleanSe[0]
                                         + res.termCleanSe[1] * res.termCleanSe[1]));
        }
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

// --- Verdicts: thresholds of the checks in the report -----------------------

enum { Good = 0, Watch = 1, Bad = 2 };

const char *levelBg[3]   = {"#d6efd6", "#fff0b8", "#f6cccc"};
const char *levelText[3] = {"OK", "CHECK", "PROBLEM"};

QString levelCell(int level){
    return QString("<td bgcolor='%1' align='center'><b>%2</b></td>").arg(levelBg[level]).arg(levelText[level]);
}

// "29.0 nm = lambda/23"
QString lamFrac(double v, double lambda){
    if (v <= 0. || lambda <= 0.)
        return nm(v) + " nm";
    return QString("%1 nm = &lambda;/%2").arg(nm(v)).arg(nm(lambda / v, 0));
}

struct Check {
    QString what;
    QString value;
    int level;
    QString hint;
};

// --- Absolute criterion: how much Strehl at 550 nm the uncertainty of this
// direction could cost.  The direction is judged on its own; averaging it with
// other rotations of the mirror can only make the error smaller.

const double kStrehlLambda = 550.;   // nm, visual
const double kStrehlOk = 0.02, kStrehlBorder = 0.05;
// rms^2 of DFTFringe's Zernike terms 4..10 over the unit disc (unnormalised polynomials)
const double kTermRms2[7] = {1. / 6, 1. / 6, 1. / 8, 1. / 8, 1. / 5, 1. / 8, 1. / 8};

// 2 sigma rms wavefront uncertainty of the low-order terms, nm.
double lowOrderUnc(const std::vector<double> &se, int first = 0, int last = 6){
    double s = 0.;
    for (int i = first; i <= last && i < (int)se.size(); ++i)
        s += se[i] * se[i] * kTermRms2[i];
    return 2. * sqrt(s);
}

double strehlLoss(double rmsNm){
    const double p = 2. * M_PI * rmsNm / kStrehlLambda;
    return 1. - exp(-p * p);
}

int strehlLevel(double loss){
    return loss <= kStrehlOk ? Good : (loss <= kStrehlBorder ? Watch : Bad);
}

struct Usability {
    int level = Good;
    double unc = 0., loss = 0.;           // all terms
    double uncAst = 0., lossAst = 0.;     // astigmatism only
    bool hasClean = false;
    int levelClean = Good;
    double uncClean = 0., lossClean = 0.;
};

Usability usability(const AccuracyResult &r){
    Usability u;
    std::vector<double> se;
    for (const AccuracyTerm &t : r.terms)
        se.push_back(t.se);
    u.unc = lowOrderUnc(se);
    u.loss = strehlLoss(u.unc);
    u.uncAst = lowOrderUnc(se, 0, 1);
    u.lossAst = strehlLoss(u.uncAst);
    u.level = strehlLevel(u.loss);
    if (r.hasClean){
        u.hasClean = true;
        u.uncClean = lowOrderUnc(r.termCleanSe);
        u.lossClean = strehlLoss(u.uncClean);
        u.levelClean = strehlLevel(u.lossClean);
    }
    return u;
}

QString usabilitySentence(const AccuracyResult &r, const Usability &u){
    const char *answer[3] = {"YES", "BORDERLINE", "NO"};
    QString s = QString("<b>Usable as one direction of the rotation test: %1.</b> ").arg(answer[u.level]);
    s += QString("The uncertainty of this series (2&sigma;, astigmatism, coma, trefoil and spherical together) is "
                 "%1 nm rms of wavefront = &lambda;/%2 at 550 nm and could cost up to %3 of Strehl "
                 "(astigmatism alone: %4 nm, %5). ")
         .arg(nm(u.unc)).arg(nm(kStrehlLambda / std::max(u.unc, 1e-6), 0)).arg(nm(u.loss, 3))
         .arg(nm(u.uncAst)).arg(nm(u.lossAst, 3));
    if (u.level == Good)
        s += "Averaging it with the other rotations only makes this smaller, so this direction will not spoil "
             "the result.";
    else if (u.level == Watch)
        s += "It can be used, but it will be one of the weaker directions of the rotation test.";
    else
        s += "This direction alone would make the rotation test unreliable.";
    if (u.hasClean && u.level != Good){
        QStringList far;
        for (const AccuracyShot &sh : r.shots)
            if (sh.flag == 2)
                far << sh.name.toHtmlEscaped();
        s += QString(" Without %1: %2 nm, %3 of Strehl - %4.")
             .arg(far.join(", ")).arg(nm(u.uncClean)).arg(nm(u.lossClean, 3))
             .arg(u.levelClean == Good ? "usable, delete the far-out shots"
                  : (u.levelClean == Watch ? "borderline" : "still not usable, take the series again"));
    } else if (u.level == Bad){
        s += " Take the series again in steadier air.";
    }
    return s;
}

std::vector<Check> qualityChecks(const AccuracyResult &r){
    std::vector<Check> c;
    const double lam = r.lambdaNm;

    // Number of shots.
    {
        Check k;
        k.what = "Number of shots";
        k.value = QString::number(r.n);
        if (r.n >= 12){
            k.level = Good;
            k.hint = "Enough shots to separate the noise from the errors that turn with the fringes.";
        } else if (r.n >= 8){
            k.level = Watch;
            k.hint = "Usable, but the noise estimate itself is rough. 12-16 shots in 6-8 fringe directions "
                     "(two per direction) give a firmer result.";
        } else {
            k.level = Bad;
            k.hint = "Fewer than 8 shots: the errors that turn with the fringes cannot be separated from the "
                     "noise. Take at least 8, better 12-16, in directions spread around the circle.";
        }
        c.push_back(k);
    }

    // Far out / deviating shots.
    {
        QStringList far;
        int dev = 0;
        for (const AccuracyShot &s : r.shots){
            if (s.flag == 2)
                far << QString("%1 (%2&deg;, %3&times;)").arg(s.name.toHtmlEscaped()).arg(nm(s.tiltDeg, 0)).arg(nm(s.ratio, 1));
            else if (s.flag == 1)
                ++dev;
        }
        Check k;
        k.what = "Shots that disagree";
        k.value = QString("%1 far out, %2 deviating").arg(far.size()).arg(dev);
        if (!far.isEmpty()){
            k.level = Bad;
            k.hint = QString("Far out: %1. Such a shot was disturbed (air, vibration, a bad fringe analysis) and "
                             "pulls the average. Delete it from the wave front list and run Average and "
                             "Accuracy again without it.").arg(far.join(", "));
            if (r.hasClean)
                k.hint += QString(" Without %1 the astigmatism would be <b>%2 &plusmn; %3 nm</b> "
                                  "instead of %4 &plusmn; %5 nm.")
                          .arg(far.size() > 1 ? "them" : "it").arg(nm(r.astigClean)).arg(nm(r.astigCleanSe)).arg(nm(r.astigMean)).arg(nm(r.astigSe));
        } else if (dev > r.n / 4){
            k.level = Watch;
            k.hint = QString("No shot is far out, but %1 of %2 deviate by more than 1.5&times; the median: "
                             "the conditions were uneven (air, vibration). The noise floor already includes them.")
                     .arg(dev).arg(r.n);
        } else {
            k.level = Good;
            k.hint = dev ? "A few shots deviate a little (up to 2.5&times; the median) - normal scatter."
                         : "All shots agree with each other.";
        }
        c.push_back(k);
    }

    // Noise floor against the RMS of the average.
    {
        const double ratio = r.rmsAvgNm > 0. ? r.floorNm / r.rmsAvgNm : 1.;
        Check k;
        k.what = "Noise floor of the average";
        k.value = QString("%1<br>%2% of the RMS (%3 nm)").arg(lamFrac(r.floorNm, lam))
                  .arg(nm(100. * ratio, 0)).arg(nm(r.rmsAvgNm));
        const bool tiny = lam > 0. && r.floorNm <= lam / 50.;
        const int need = r.rmsAvgNm > 0. ? int(ceil(r.n * pow(3. * r.floorNm / r.rmsAvgNm, 2))) : 0;
        if (ratio <= 1. / 3. || tiny){
            k.level = Good;
            k.hint = tiny && ratio > 1. / 3.
                   ? "The noise floor is very small in absolute terms (below &lambda;/50); the mirror itself is very smooth."
                   : QString("The average stands %1&times; above its own noise: the shape you see is real, "
                             "it is not air or vibration.").arg(nm(1. / std::max(ratio, 1e-6), 1));
        } else if (ratio <= 0.5 || (lam > 0. && r.floorNm <= lam / 30.)){
            k.level = Watch;
            k.hint = QString("Part of the fine detail of the average is still noise. About %1 shots in total "
                             "would bring the floor down to a third of the RMS (the floor falls with the "
                             "square root of the number of shots).").arg(need);
        } else {
            k.level = Bad;
            k.hint = QString("The average is largely noise. Take more shots (about %1 in total) or measure "
                             "in steadier air.").arg(need);
        }
        c.push_back(k);
    }

    // Independence of the shots.
    {
        const double q = r.n > 0 ? r.nEff / r.n : 0.;
        Check k;
        k.what = "Independent shots";
        k.value = QString("%1 of %2 (correlation %3)").arg(nm(r.nEff)).arg(r.n).arg(nm(r.corr, 2));
        if (q >= 0.85){
            k.level = Good;
            k.hint = "Each shot saw different air - every shot counts fully.";
        } else if (q >= 0.7){
            k.level = Watch;
            k.hint = "Consecutive shots are somewhat alike (slow air drift). The noise floor already takes it "
                     "into account; leave more time between shots.";
        } else {
            k.level = Bad;
            k.hint = "Consecutive shots repeat the same error (slow air current or drift): more shots in a row "
                     "help little. Let the air settle, leave more time between shots, or block the air flow.";
        }
        c.push_back(k);
    }

    // Fringe directions.
    {
        Check k;
        k.what = "Fringe directions";
        const double bal = std::max(r.balance1, r.balance2);
        k.value = QString("gap %1&deg;, balance %2").arg(nm(r.maxGapDeg, 0)).arg(nm(bal, 2));
        if (r.fringeModel)
            k.value += QString("<br>%1 nm per shot, %2 nm left").arg(nm(r.lockedNm)).arg(nm(r.biasNm));
        if (r.maxGapDeg > 120.){
            k.level = Bad;
            k.hint = "The directions leave a large gap: the error that turns with the fringes does not cancel "
                     "and may add astigmatism to the average that the noise figures do not show. Spread the "
                     "fringe directions around the whole circle (e.g. every 45&deg;).";
        } else if (!r.fringeModel){
            k.level = Watch;
            k.hint = "The error that turns with the fringes could not be estimated (too few shots or too little "
                     "tilt), so it is unknown how much of it is left in the average.";
        } else if (r.biasNm > r.floorNm){
            k.level = Bad;
            k.hint = QString("The directions are unbalanced: %1 nm of the turning error stays in the average, more "
                             "than the noise floor. Add shots in the missing directions.").arg(nm(r.biasNm));
        } else if (r.maxGapDeg > 60. || bal > 0.15 || r.biasNm > r.floorNm / 3.){
            k.level = Watch;
            k.hint = QString("The turning error mostly cancels (%1 nm left), but the directions are not spread "
                             "evenly. Shots every 45&deg; (or 60&deg;) around the circle are best.").arg(nm(r.biasNm));
        } else {
            k.level = Good;
            k.hint = QString("Each shot carries %1 nm of error that turns with the fringes, but the directions are "
                             "spread evenly and only %2 nm of it is left in the average.")
                     .arg(nm(r.lockedNm)).arg(nm(r.biasNm));
        }
        c.push_back(k);
    }

    // Stability of the astigmatism (leave one out).
    {
        const double span = r.looAstigMax - r.looAstigMin;
        const double ref = std::max(r.astigMean, lam > 0. ? lam / 20. : 1.);
        const double rel = span / ref;
        QString who = (r.looAstigShot >= 0 && r.looAstigShot < (int)r.shots.size())
                    ? r.shots[r.looAstigShot].name.toHtmlEscaped() : QString();
        Check k;
        k.what = "Astigmatism, leaving out one shot";
        k.value = QString("%1 &ndash; %2 nm").arg(nm(r.looAstigMin)).arg(nm(r.looAstigMax));
        if (rel <= 0.15){
            k.level = Good;
            k.hint = "No single shot decides the astigmatism.";
        } else if (rel <= 0.3){
            k.level = Watch;
            k.hint = QString("One shot moves the astigmatism noticeably (most: %1). Look at it in the Shots "
                             "table and in the astigmatism plot.").arg(who);
        } else {
            k.level = Bad;
            k.hint = QString("The astigmatism depends strongly on one shot (%1). If it stands out in the Shots "
                             "table or the plot, delete it from the analysis; otherwise take more shots before "
                             "trusting the astigmatism.").arg(who);
        }
        c.push_back(k);
    }
    return c;
}

// Significance of a Zernike term: how many standard errors from zero.
QString significanceCell(double mean, double se){
    if (se <= 0.)
        return "<td></td>";
    const double z = fabs(mean) / se;
    if (z >= 3.)
        return QString("<td bgcolor='#d6efd6'>real (%1&sigma;)</td>").arg(nm(z, 1));
    if (z >= 2.)
        return QString("<td bgcolor='#fff0b8'>probably real (%1&sigma;)</td>").arg(nm(z, 1));
    return QString("<td>noise level (%1&sigma;)</td>").arg(nm(z, 1));
}

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

    // Verdict.
    {
        const std::vector<Check> checks = qualityChecks(r);
        const Usability use = usability(r);
        int worst = use.level;
        for (const Check &k : checks)
            worst = std::max(worst, k.level);
        const char *overall[3] = {
            "GOOD - the average can be trusted to within its noise floor.",
            "USABLE - look at the items marked CHECK.",
            "NOT RELIABLE YET - fix the items marked PROBLEM (e.g. delete far-out shots) before trusting the average."};
        h += "<h3>Verdict</h3>";
        h += QString("<table border='1' cellspacing='0' cellpadding='6' width='100%'><tr><td bgcolor='%1'>"
                     "%2</td></tr></table><p></p>").arg(levelBg[use.level]).arg(usabilitySentence(r, use));
        h += QString("<table border='1' cellspacing='0' cellpadding='6' width='100%'><tr><td bgcolor='%1'>"
                     "<b>%2</b></td></tr></table>").arg(levelBg[worst]).arg(overall[worst]);
        h += "<p></p><table border='1' cellspacing='0' cellpadding='4' width='100%'>";
        h += "<tr><th>Check</th><th>Value</th><th></th><th>What it means / what to do</th></tr>";
        for (const Check &k : checks)
            h += QString("<tr><td>%1</td><td>%2</td>%3<td>%4</td></tr>")
                 .arg(k.what).arg(k.value).arg(levelCell(k.level)).arg(k.hint);
        h += "</table>";
        h += "<p><i>Not covered by these checks:</i> errors that are the same in every shot - the stand, the "
             "reference, the null, steady air currents. A term can be &quot;real&quot; here and still come from "
             "the test setup. To tell: rotate the mirror on the stand (90&deg; or more) and measure again - "
             "astigmatism or coma that turns with the mirror belongs to the mirror, what stays in place belongs "
             "to the setup.</p>";
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
    h += "<tr><th>Term</th><th>Average</th><th>&plusmn; std. error</th><th>single shot scatter</th><th>Is it real?</th>";
    if (r.hasClean)
        h += QString("<th>without far-out shots (%1)</th>").arg(r.nClean);
    h += "</tr>";
    for (size_t i = 0; i < r.terms.size(); ++i){
        const AccuracyTerm &t = r.terms[i];
        h += QString("<tr><td>%1</td><td align='right'>%2</td><td align='right'>%3</td><td align='right'>%4</td>%5")
             .arg(t.name).arg(nm(t.mean)).arg(nm(t.se)).arg(nm(t.sigma1)).arg(significanceCell(t.mean, t.se));
        if (r.hasClean && i < r.termCleanMean.size())
            h += QString("<td align='right'>%1 &plusmn; %2</td>").arg(nm(r.termCleanMean[i])).arg(nm(r.termCleanSe[i]));
        h += "</tr>";
    }
    h += "</table>";
    h += "<p>Spherical is before the null. The std. error is 1 sigma; about 95% of repeated averages fall within "
         "2 std. errors. <b>Is it real?</b> compares the term with its std. error: 3&sigma; or more - real; "
         "2-3&sigma; - probably real; below 2&sigma; - it cannot be told from noise, treat it as zero. "
         "If the single shot scatter of a term is much larger than in the other terms, look for a far-out shot.</p>";
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
         "Deviating: over 1.5&times; the median. Far out: 2.5&times; or more - delete it from the analysis. "
         "No shot is left out of the average automatically.</p>";

    // Reading guide.
    h += "<h3>How to read this report</h3>";
    h += "<table border='1' cellspacing='0' cellpadding='4' width='100%'>";
    h += "<tr><th>Item</th><th>What it is</th><th>Good</th><th>Not good</th></tr>";
    auto guide = [&](const char *item, const char *what, const char *good, const char *bad){
        h += QString("<tr><td><b>%1</b></td><td>%2</td><td bgcolor='#eef8ee'>%3</td><td bgcolor='#fbecec'>%4</td></tr>")
             .arg(item).arg(what).arg(good).arg(bad);
    };
    guide("Usable as one direction of the rotation test",
          "The std. errors of astigmatism, coma, trefoil and spherical, turned into one rms wavefront figure "
          "(2&sigma;) and into the Strehl ratio at 550 nm it could cost (Mar&eacute;chal). The series is judged "
          "on its own, without relying on the stand or on the other rotations: averaging over rotations can only "
          "make the error smaller.",
          "YES: up to 0.02 of Strehl (about 12 nm rms, &lambda;/45 at 550 nm). For astigmatism alone this is a "
          "std. error up to about 11 nm per term.",
          "BORDERLINE: 0.02-0.05 (12-20 nm rms) - usable, but a weak direction. NO: above 0.05 - delete far-out "
          "shots (the report says if that is enough) or take the series again.");
    guide("Attention",
          "Problems found automatically: far-out shots, gaps in the fringe directions, correlated shots.",
          "Empty.",
          "Anything listed - deal with it first, the rest of the report is less reliable until then.");
    guide("Noise floor of the average",
          "How far the average itself may still be from the truth because of random errors (air, vibration, "
          "fringe analysis). It is the scatter of one shot divided by the square root of the number of "
          "independent shots: 4&times; more shots halve it.",
          "Below a third of the RMS of the average, or below &lambda;/50. The shape of the average is then real.",
          "Above half of the RMS: the fine detail of the average is mostly noise. Take more shots or wait for "
          "steadier air.");
    guide("Scatter of a single shot",
          "How much one shot differs from the next. The random part is what is left after the error that turns "
          "with the fringes; about the plain average it is a little larger.",
          "Small compared with the RMS of the average. The two scatter figures close to each other.",
          "Much larger than the RMS: unsteady air or vibration. The plain one much larger than the random one: "
          "a strong error turning with the fringes.");
    guide("Independent shots",
          "Shots taken one after another may see the same air and repeat the same error. The correlation of "
          "consecutive shots turns into an effective number of independent shots.",
          "Correlation around 0 or negative (up to about +0.08); 85% or more of the shots independent.",
          "Correlation above about +0.18; fewer than 70% of the shots independent. Slow air current: more shots "
          "in a row do not help, time between them does.");
    guide("Error turning with the fringes",
          "Part of the error of each shot follows the direction of the fringes (mostly astigmatism from the "
          "interferometer and the analysis). In shots spread evenly around the circle it cancels.",
          "Whatever it is per shot, the part left in the average is well below the noise floor.",
          "The part left in the average is comparable with or larger than the noise floor: add shots in the "
          "missing directions.");
    guide("Balance and largest gap",
          "How evenly the fringe directions cover the circle. Balance 0 means perfectly even.",
          "Balance below 0.15, largest gap 60&deg; or less (e.g. 8 directions every 45&deg;, two shots in each).",
          "Gap over 120&deg; - the turning error cannot even be estimated. Balance above 0.3.");
    guide("RMS of the average",
          "The surface error of the average without piston, tilt, focus and spherical (the null is not applied). "
          "Only a reference for the noise floor - the full analysis is in the main report.",
          "-", "-");
    guide("Zernike terms",
          "Each term of the average with its std. error (1&sigma;) and the scatter of one shot.",
          "Average at least 3 std. errors from zero - real. Small std. error compared with the tolerance you "
          "work to.",
          "Below 2 std. errors - indistinguishable from zero. Std. error comparable with the tolerance - the "
          "term is not measured well enough.");
    guide("Astigmatism, leaving out one shot",
          "The astigmatism of the average recomputed once without each shot. Shows whether one shot "
          "decides the result.",
          "The range is narrow - within about 15% of the astigmatism.",
          "A wide range, especially if one end is far from the rest: one shot dominates. Find it in the Shots "
          "table (its Astig X/Y stands out) and in the plot (a mark far from the cloud).");
    guide("Astigmatism plot",
          "Each mark is one shot (Astig X against Y), the cross is the average, the dashed circle is 2 std. errors.",
          "The marks form one cloud around the cross; blue marks.",
          "A mark far outside the cloud (orange - deviating, red - far out).");
    guide("Shots: deviation and &times; median",
          "What is left of each shot after the average and the turning error, compared with the typical shot.",
          "Up to 1.5&times; the median.",
          "1.5-2.5&times; - deviating, tolerable. 2.5&times; or more - far out: delete it from the wave front "
          "list and run Average and Accuracy again. Check afterwards that the fringe directions are still "
          "balanced.");
    h += "</table>";
    h += "<p>Units: all values are nm of <b>wavefront</b>. The surface error is half of it. "
         "&lambda;/x figures use the test wavelength.</p>";
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
    // Not modal and minimizable: it can stay open next to the wave front list.
    QDialog *pdlg = new QDialog(parent, Qt::Window | Qt::WindowTitleHint | Qt::WindowSystemMenuHint |
                                Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint |
                                Qt::WindowCloseButtonHint);
    pdlg->setAttribute(Qt::WA_DeleteOnClose);
    QDialog &dlg = *pdlg;
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
    dlg.show();
    return written;
}
