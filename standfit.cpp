/******************************************************************************
**
**  standfit.cpp -- see standfit.h for the model and the closed form.
**
******************************************************************************/
#include "standfit.h"

#include <algorithm>
#include <cmath>

namespace {

const double PI = 3.14159265358979323846;

// A rotation of the pair (x,y) by phi.
inline void rot(double phi, double x, double y, double &ox, double &oy)
{
    const double c = cos(phi), s = sin(phi);
    ox = c * x - s * y;
    oy = s * x + c * y;
}

// Median of a copy of v.  v is small (one entry per rotation).
double median(std::vector<double> v)
{
    if (v.empty())
        return 0.;
    std::sort(v.begin(), v.end());
    const size_t h = v.size() / 2;
    return (v.size() % 2) ? v[h] : 0.5 * (v[h - 1] + v[h]);
}

} // namespace

void zernikeOrderPairs(int nTerms, std::vector<StandFitOrder> &out)
{
    out.clear();
    // Wyant ordering is grouped in blocks: block b covers the indices
    // [b*b, b*b + 2*b].  Inside it the pairs run from m = b down to m = 1,
    // two indices each, and the last index of the block is the m = 0 term.
    //     b=2: 4,5 astig (m=2)   6,7 coma (m=1)   8 spherical (m=0)
    //     b=3: 9,10 trefoil      11,12 2nd astig  13,14 2nd coma  15 m=0
    for (int b = 1; b * b < nTerms; ++b) {
        for (int t = 0; t < b; ++t) {
            StandFitOrder o;
            o.m      = b - t;
            o.cosNdx = b * b + 2 * t;
            o.sinNdx = o.cosNdx + 1;
            if (o.sinNdx >= nTerms)
                continue;
            if (o.cosNdx == 1)
                continue;           // tilt: alignment, not stand deformation
            out.push_back(o);
        }
    }
}

StandFit fitStandZernikes(const std::vector<double> &anglesDeg,
                          const std::vector<std::vector<double> > &zerns,
                          int nTerms)
{
    StandFit fit;
    const int N = static_cast<int>(anglesDeg.size());
    fit.n = N;
    fit.standZerns.assign(nTerms, 0.);
    fit.mirrorSE.assign(nTerms, 0.);
    fit.residual.assign(N > 0 ? N : 0, 0.);
    fit.residualAstig.assign(N > 0 ? N : 0, 0.);
    fit.outlier.assign(N > 0 ? N : 0, false);

    if (N < 2 || static_cast<int>(zerns.size()) != N)
        return fit;                      // one rotation says nothing

    zernikeOrderPairs(nTerms, fit.orders);
    if (fit.orders.empty())
        return fit;

    // Residual accumulators.  sumSq is over every component actually fitted,
    // which is what sigma and the per rotation residuals are built from.
    std::vector<double> resSq(N, 0.);
    std::vector<double> resSqLow(N, 0.);
    std::vector<double> resSqAstig(N, 0.);
    std::vector<double> sumSqOrder(fit.orders.size(), 0.);
    double sumSq = 0.;
    int    nFittedPairs = 0;
    int    nLowPairs = 0;       // primary astig, coma and trefoil

    for (size_t k = 0; k < fit.orders.size(); ++k) {
        StandFitOrder &o = fit.orders[k];
        const int ci = o.cosNdx, si = o.sinNdx;

        double cr = 0., ci_ = 0.;        // c = sum exp(i*m*theta)
        double ux = 0., uy = 0.;         // u = sum a_i
        double vx = 0., vy = 0.;         // v = sum R(m*theta)^T a_i
        for (int i = 0; i < N; ++i) {
            const double phi = o.m * anglesDeg[i] * PI / 180.;
            const double ax = (ci < static_cast<int>(zerns[i].size())) ? zerns[i][ci] : 0.;
            const double ay = (si < static_cast<int>(zerns[i].size())) ? zerns[i][si] : 0.;
            cr += cos(phi);
            ci_ += sin(phi);
            ux += ax;
            uy += ay;
            double tx, ty;
            rot(-phi, ax, ay, tx, ty);   // R^T = R(-phi)
            vx += tx;
            vy += ty;
        }

        const double cmag2 = cr * cr + ci_ * ci_;
        const double den   = static_cast<double>(N) * N - cmag2;
        o.balance = (N > 0) ? sqrt(cmag2) / N : 1.;

        // den/N^2 is how much of the information survives.  Below a few
        // percent the two unknowns are numerically the same thing.
        if (den <= 0.05 * static_cast<double>(N) * N) {
            o.separable = false;
            continue;
        }
        o.separable = true;
        ++nFittedPairs;
        const bool lowOrder = (ci == 4 || ci == 6 || ci == 9);
        if (lowOrder)
            ++nLowPairs;

        // M = (N*v - C^T*u) / den,  C^T*u = R(-arg c)*|c|*u written out
        double ctux, ctuy;
        ctux = cr * ux + ci_ * uy;       // C^T = [[cr, ci],[-ci, cr]]
        ctuy = -ci_ * ux + cr * uy;
        o.mirrorX = (N * vx - ctux) / den;
        o.mirrorY = (N * vy - ctuy) / den;

        // S = (u - C*M) / N
        const double cmx = cr * o.mirrorX - ci_ * o.mirrorY;
        const double cmy = ci_ * o.mirrorX + cr * o.mirrorY;
        o.standX = (ux - cmx) / N;
        o.standY = (uy - cmy) / N;

        fit.standZerns[ci] = o.standX;
        fit.standZerns[si] = o.standY;

        for (int i = 0; i < N; ++i) {
            const double phi = o.m * anglesDeg[i] * PI / 180.;
            const double ax = (ci < static_cast<int>(zerns[i].size())) ? zerns[i][ci] : 0.;
            const double ay = (si < static_cast<int>(zerns[i].size())) ? zerns[i][si] : 0.;
            double rx, ry;
            rot(phi, o.mirrorX, o.mirrorY, rx, ry);
            const double ex = ax - o.standX - rx;
            const double ey = ay - o.standY - ry;
            const double e2 = ex * ex + ey * ey;
            resSq[i] += e2;
            if (lowOrder)
                resSqLow[i] += e2;
            if (o.m == 2 && ci == 4)
                resSqAstig[i] = e2;
            sumSqOrder[k] += e2;
            sumSq += e2;
        }
    }

    if (nFittedPairs == 0)
        return fit;

    // Each fitted pair costs 4 parameters and brings 2*N measurements.
    fit.dof = 2 * N * nFittedPairs - 4 * nFittedPairs;
    fit.sigma = (fit.dof > 0) ? sqrt(sumSq / fit.dof) : 0.;
    fit.valid = true;

    // The residual that is reported and that outliers are judged by covers the
    // primary astigmatism, coma and trefoil only.  Those are the terms a test
    // stand actually bends into the glass; averaging in twenty high orders
    // would bury a rotation that went wrong in terms that carry no stand
    // signal at all.
    for (int i = 0; i < N; ++i) {
        fit.residual[i] = (nLowPairs > 0)
                ? sqrt(resSqLow[i] / (2. * nLowPairs))
                : sqrt(resSq[i] / (2. * nFittedPairs));
        fit.residualAstig[i] = sqrt(resSqAstig[i]);
    }

    // Standard error of the mirror terms.  The normal matrix block for M
    // inverts to N/(N^2-|c|^2) * I, so the error depends only on the angle set
    // and on the scatter.  The scatter is taken PER ORDER: low orders carry
    // most of the seeing and of the real stand movement, so one common sigma
    // would overstate the error of the high orders and understate the astig.
    const int dofOrder = 2 * N - 4;
    for (size_t k = 0; k < fit.orders.size(); ++k) {
        const StandFitOrder &o = fit.orders[k];
        if (!o.separable)
            continue;
        const double cmag2 = o.balance * o.balance * static_cast<double>(N) * N;
        const double den   = static_cast<double>(N) * N - cmag2;
        const double sigmaK = (dofOrder > 0) ? sqrt(sumSqOrder[k] / dofOrder) : 0.;
        const double se     = sigmaK * sqrt(static_cast<double>(N) / den);
        fit.mirrorSE[o.cosNdx] = se;
        fit.mirrorSE[o.sinNdx] = se;
    }

    // Outliers.  With four rotations the fit still has something to say; with
    // fewer, dropping one would leave it underdetermined, so nothing is
    // flagged.  Robust z score on the per rotation residual (MAD based) keeps
    // one bad rotation from setting the scale it is judged against.
    if (N >= 5) {
        const double med = median(fit.residual);
        std::vector<double> dev(N);
        for (int i = 0; i < N; ++i)
            dev[i] = fabs(fit.residual[i] - med);
        const double mad = median(dev);
        if (mad > 1e-9) {
            for (int i = 0; i < N; ++i) {
                const double z = 0.6745 * (fit.residual[i] - med) / mad;
                if (z > 3.5)
                    fit.outlier[i] = true;
            }
        }
    }

    return fit;
}
