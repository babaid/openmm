#include <string.h>
#include <sstream>
#include<cmath>
#include "SimTKOpenMMUtilities.h"
#include "ReferenceCutoffAngleBondIxn.h"
#include "ReferenceForce.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceAngleBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCutoffAngleBondIxn::ReferenceCutoffAngleBondIxn() : usePeriodic(false) {
}

/**---------------------------------------------------------------------------------------

   ReferenceAngleBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCutoffAngleBondIxn::~ReferenceCutoffAngleBondIxn() {
}

void ReferenceCutoffAngleBondIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Get dEdR and energy term for angle bond

   @param  cosine               cosine of angle
   @param  angleParameters      angleParameters: angleParameters[0] = angle in radians
                                                 angleParameters[1] = k (force constant)
   @param  dEdR                 output dEdR
   @param  energyTerm           output energyTerm

   --------------------------------------------------------------------------------------- */

void ReferenceCutoffAngleBondIxn::getPrefactorsGivenAngleCosine(double cosine, vector<double>& angleParameters,
                                                          double* dEdR, double* energyTerm) const {

    double angle;
    if (cosine >= 1.0) {
        angle = 0.0;
    } else if (cosine <= -1.0) {
        angle = PI_M;
    } else {
        angle = acos(cosine);
    }
    double deltaIdeal         = angle - angleParameters[0];
    double deltaIdeal2        = deltaIdeal*deltaIdeal;

    *dEdR                          = angleParameters[1]*deltaIdeal;
    *energyTerm                    = 0.5*angleParameters[1]*deltaIdeal2;

}

/**---------------------------------------------------------------------------------------

   Calculate Angle Bond ixn

   @param atomIndices      two bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters: parameters[0] = ideal bond length
                                       parameters[1] = bond k (includes factor of 2)
   @param forces           force array (forces added)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCutoffAngleBondIxn::calculateBondIxn(vector<int>& atomIndices,
                                             vector<Vec3>& atomCoordinates,
                                             vector<double>& parameters,
                                             vector<Vec3>& forces,
                                             double* totalEnergy, double* energyParamDerivs) {

    static const int LastAtomIndex      = 3;

    double deltaR[2][ReferenceForce::LastDeltaRIndex];

    // ---------------------------------------------------------------------------------------

    // get deltaR, R2, and R between 2 atoms
    int atomAIndex = atomIndices[0];
    int atomBIndex = atomIndices[1];
    int atomCIndex = atomIndices[2];
    if (usePeriodic) {
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], boxVectors, deltaR[0]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomCIndex], atomCoordinates[atomBIndex], boxVectors, deltaR[1]);
    }
    else {
        ReferenceForce::getDeltaR(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR[0]);
        ReferenceForce::getDeltaR(atomCoordinates[atomCIndex], atomCoordinates[atomBIndex], deltaR[1]);
    }

    double pVector[3];
    SimTKOpenMMUtilities::crossProductVector3(deltaR[0], deltaR[1], pVector);
    double rp = sqrt(DOT3(pVector, pVector));
    if (rp < 1.0e-06) {
        rp = 1.0e-06;
    }
    double dot             = DOT3(deltaR[0], deltaR[1]);
    double cosine          = dot/sqrt((deltaR[0][ReferenceForce::R2Index]*deltaR[1][ReferenceForce::R2Index]));

    double dEdR;
    double energy;
    getPrefactorsGivenAngleCosine(cosine, parameters, &dEdR, &energy);
    dEdR*=parameters[2];
    energy*=parameters[2];
    double termA           =  dEdR/(deltaR[0][ReferenceForce::R2Index]*rp);
    double termC           = -dEdR/(deltaR[1][ReferenceForce::R2Index]*rp);

    double deltaCrossP[LastAtomIndex][3];
    SimTKOpenMMUtilities::crossProductVector3(deltaR[0], pVector, deltaCrossP[0]);
    SimTKOpenMMUtilities::crossProductVector3(deltaR[1], pVector, deltaCrossP[2]);

    for (int ii = 0; ii < 3; ii++) {
        deltaCrossP[0][ii] *= termA;
        deltaCrossP[2][ii] *= termC;
        deltaCrossP[1][ii]  = -(deltaCrossP[0][ii] + deltaCrossP[2][ii]);
    }

    // accumulate forces

    for (int jj = 0; jj < LastAtomIndex; jj++) {
        for (int ii = 0; ii < 3; ii++) {
            forces[atomIndices[jj]][ii] += deltaCrossP[jj][ii];
        }
    }

    // accumulate energies

    if (totalEnergy != NULL)
        *totalEnergy += energy;
}
