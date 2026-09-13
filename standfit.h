/******************************************************************************
**
**  standfit.h  --  least squares separation of test stand and mirror terms
**                  from a set of measurements taken at different rotations.
**
**  Part of the DFTFringe-AR fork.  Licensed like the rest of DFTFringe
**  (GNU GPL v3).
**
**  The model
**  ---------
**  A measurement taken with the mirror rotated by theta is
**
**      a_i = S + R(m*theta_i) * M
**
**  where a_i is the (X,Y) pair of Zernike coefficients of one azimuthal
**  order m, S is the term the test stand contributes (fixed in the lab
**  frame) and M is the term belonging to the mirror (fixed in the mirror
**  frame, so it rotates with it).  R is the usual 2x2 rotation.
**
**  Averaging counter rotated wavefronts - what DFTFringe did before - is
**  the special case of this fit for an angle set that happens to be
**  balanced, i.e. sum(exp(i*m*theta)) = 0.  The pair 0/90 deg balances
**  order m=2 (astigmatism) only: stand coma (m=1) and trefoil (m=3) pass
**  straight through into the "stand removed" result.  Solving the model
**  instead of averaging removes the stand for every separable order and,
**  more importantly, leaves residuals that say how badly the assumption
**  "the stand does the same thing at every rotation" is violated.
**
**  Closed form
**  -----------
**  With c = sum(exp(i*m*theta_i)), C = [[Re c, -Im c],[Im c, Re c]],
**  u = sum(a_i) and v = sum(R(m*theta_i)^T a_i), the normal equations are
**
**      N*S + C*M = u
**      C^T*S + N*M = v
**
**  and because C^T*C = |c|^2 * I they collapse to
**
**      M = (N*v - C^T*u) / (N^2 - |c|^2)
**      S = (u - C*M) / N
**
**  The denominator is the separability of the angle set for that order:
**  N^2 - |c|^2 = 0 means every rotation looks the same to order m and
**  nothing can be told apart.  m = 0 (defocus, spherical) is that case
**  for EVERY angle set - rotation can never separate a rotationally
**  symmetric stand deformation from the mirror figure.
**
******************************************************************************/
#ifndef STANDFIT_H
#define STANDFIT_H

#include <vector>

struct StandFitOrder {
    int    m = 0;             // azimuthal order
    int    cosNdx = 0;        // Wyant index of the X (cos) term
    int    sinNdx = 0;        // Wyant index of the Y (sin) term
    double standX = 0.;       // stand contribution, lab frame
    double standY = 0.;
    double mirrorX = 0.;      // mirror contribution, mirror frame
    double mirrorY = 0.;
    double balance = 1.;      // |sum exp(i*m*theta)| / N.  0 = ideal, 1 = blind
    bool   separable = false; // false -> nothing subtracted for this order
};

struct StandFit {
    int    n = 0;                        // number of rotations
    int    dof = 0;                      // degrees of freedom of the fit
    double sigma = 0.;                   // rms residual per component, waves
    bool   valid = false;
    std::vector<StandFitOrder> orders;   // separable and non separable alike
    std::vector<double> standZerns;      // size nTerms, only separable m>0 filled
    std::vector<double> mirrorSE;        // size nTerms, standard error of M
    std::vector<double> residual;        // per rotation, rms over primary astig/coma/trefoil
    std::vector<double> residualAstig;   // per rotation, astig only (m=2)
    std::vector<bool>   outlier;         // per rotation, left out of the average
    std::vector<bool>   suspect;         // per rotation, worth measuring again
};

// anglesDeg[i]  - rotation of the mirror for measurement i, degrees
// zerns[i]      - that measurement's Zernike coefficients (Wyant order)
// nTerms        - how many coefficients to consider (Z_TERMS)
StandFit fitStandZernikes(const std::vector<double> &anglesDeg,
                          const std::vector<std::vector<double> > &zerns,
                          int nTerms);

// Which Wyant index pairs belong to which azimuthal order.  Piston, tilt
// and every m = 0 term are left out: tilt is alignment rather than stand
// deformation, and m = 0 cannot be separated by rotation at all.
void zernikeOrderPairs(int nTerms, std::vector<StandFitOrder> &out);

#endif // STANDFIT_H
