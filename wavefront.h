/******************************************************************************
**
**  Copyright 2016 Dale Eason
**  This file is part of DFTFringe
**  is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation version 3 of the License

** DFTFringe is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with DFTFringe.  If not, see <http://www.gnu.org/licenses/>.

****************************************************************************/
#ifndef WAVEFRONT_H
#define WAVEFRONT_H
#include <opencv2/opencv.hpp>
#include "Circleoutline.h"
#include <QPointF>
#include <QString>

enum class WavefrontOrigin : std::uint8_t {
    Unknown = 0,
    Igram,
    File,
    Simulation,
    Demo,
    Average,
    Subtraction,
    Zernikes,
    Smoothed
};
// How much of a rotation average (stand removal) is noise.  Worked out by
// SurfaceManager::computeRotationNoise() from the rotations that went into
// the average, saved in the .wft footer so a report made later from the
// saved file can still show it.  Lengths are nm: they do not depend on the
// output wavelength the report happens to use.
struct RotationNoise {
    bool   valid = false;
    int    n = 0;                    // rotations in the average
    double sigmaNm = 0.;             // rms noise of a single rotation
    double floorNm = 0.;             // noise left in the average, sigma/sqrt(n)
    double rmsNm = 0.;               // rms of the average when this was worked out
    QString removed;                 // terms the metric left out then, "0,1,2,3,8"
    QString dropped;                 // rotations the stand fit left out, "0,180"
    std::vector<double> angles;      // per rotation in the average
    std::vector<double> residualNm;
    std::vector<bool>   deviating;
    // The same with the surface smoothed afterwards: sigma of one rotation
    // and rms of the average for each Gaussian blur percent ...
    std::vector<double> blurPct, blurSigmaNm, blurRmsNm;
    // ... and for each order of Zernike smoothing.
    std::vector<int>    zernOrder;
    std::vector<double> zernSigmaNm, zernRmsNm;
    int smoothOrder = 0;             // > 0: this wavefront is a Zernike smoothing to that order
};

class wavefront { public:
    wavefront();
    ~wavefront() = default;
    wavefront( const wavefront &wf);
    wavefront( wavefront && );
    wavefront& operator=( const wavefront & ) = default;
    wavefront& operator=(wavefront &&) = default;
    void cloneMatricesFrom(const wavefront &wf);

    cv::Mat_<double> data;
    cv::Mat_<double> nulledData;
    cv::Mat_<uint8_t> mask;
    cv::Mat_<double> workData;
    cv::Mat_<uint8_t> workMask;
    std::vector<double> InputZerns;

    QVector<std::vector<cv::Point> > regions;

    std::vector<bool> zernEnablesApplied;

    double gaussian_diameter;
                  CircleOutline m_outside;
                  CircleOutline m_inside;
    double GBSmoothingValue;
    double gbValue;
    WavefrontOrigin m_origin;
    QString name;
    double lambda;
    double diameter;
    double roc;
    double min;
    double max;
    double std;
    double mean;
    bool dirtyZerns;
    bool regions_have_been_expanded;
    bool gbEnabled;
    bool wasSmoothed;
    bool useSANull;
    bool m_manuallyInverted; // true if user inverted this wavefront
    RotationNoise rotationNoise;
    QString sourceFolder;   // folder of the interferogram or .wft it came from (not saved)
};

#endif // WAVEFRONT_H
