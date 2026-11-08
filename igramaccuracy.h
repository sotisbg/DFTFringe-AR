#ifndef IGRAMACCURACY_H
#define IGRAMACCURACY_H

// How reliable is the average of several interferograms of one mirror taken in
// one direction (the mirror is not rotated, the fringes are turned between
// shots)?  See igramaccuracy.cpp for the model.

#include <QList>
#include <QString>
#include <QStringList>
#include <vector>

class wavefront;
class QWidget;

struct AccuracyShot {
    QString name;
    double tiltDeg = 0.;      // direction of the tilt (fringe orientation), degrees
    double tiltWaves = 0.;    // amount of tilt, waves
    double residualNm = 0.;   // RMS of this shot after the model, nm wavefront
    double ratio = 0.;        // residual / median residual
    int flag = 0;             // 0 ok, 1 deviating, 2 far out
    double z[11] = {0};       // Zernike terms 0..10 of this shot, nm wavefront
};

struct AccuracyTerm {
    QString name;
    double mean = 0.;         // what Average gives, nm
    double se = 0.;           // its standard error, nm
    double sigma1 = 0.;       // scatter of a single shot, nm
};

struct AccuracyResult {
    bool valid = false;
    QString error;
    QStringList warnings;
    int n = 0;
    double lambdaNm = 0.;
    bool fringeModel = false;     // the terms that turn with the fringes were fitted
    double sigmaPlainNm = 0.;     // scatter of a single shot about the plain average
    double sigmaNm = 0.;          // scatter left after the fringe model
    double lockedNm = 0.;         // RMS of the part that turns with the fringes, per shot
    double corr = 0.;             // correlation of consecutive residuals
    double nEff = 0.;             // effective number of independent shots
    double floorNm = 0.;          // noise floor of the average
    double biasNm = 0.;           // shift of the plain average from unbalanced angles
    double balance1 = 0., balance2 = 0.;   // |mean e^(i k phi)|, k = 1, 2
    double maxGapDeg = 0.;        // largest gap between tilt directions
    double rmsAvgNm = 0.;         // RMS of the average, piston/tilt/focus/spherical removed
    std::vector<AccuracyShot> shots;
    std::vector<AccuracyTerm> terms;
    double astigMean = 0., astigCorrected = 0., astigSe = 0.;   // |astigmatism|
    double astigLockedNm = 0.;    // astigmatism turning with the fringes
    double looAstigMin = 0., looAstigMax = 0.;
    double looSphMin = 0., looSphMax = 0.;
    int looAstigShot = -1;        // the shot whose removal moves the astigmatism most
    double astigAngleDeg = 0.;
    // The same without the far-out shots (what Average would give without them).
    bool hasClean = false;
    int nClean = 0;
    double astigClean = 0., astigCleanSe = 0.;
    std::vector<double> termCleanMean, termCleanSe;   // same order as terms
};

AccuracyResult computeIgramAccuracy(const QList<wavefront *> &wfs);

// HTML of the report; the astigmatism plot is referenced as "accuracy://astig".
QString accuracyReportHtml(const AccuracyResult &r, const QString &title);

// Writes the report as a PDF and shows it.  Returns the PDF path ("" if not written).
QString showAccuracyReport(const AccuracyResult &r, const QString &title,
                           const QString &pdfPath, QWidget *parent);

#endif // IGRAMACCURACY_H
