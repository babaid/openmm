#include <string.h>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceCutoffProperDihedralBond.h"
#include "ReferenceForce.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceProperDihedralBond constructor

   --------------------------------------------------------------------------------------- */

ReferenceCutoffProperDihedralBond::ReferenceCutoffProperDihedralBond() : usePeriodic(false) {
}

/**---------------------------------------------------------------------------------------

   ReferenceProperDihedralBond destructor

   --------------------------------------------------------------------------------------- */

ReferenceCutoffProperDihedralBond::~ReferenceCutoffProperDihedralBond() {
}

void ReferenceCutoffProperDihedralBond::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Calculate proper dihedral bond ixn

   @param atomIndices      atom indices of 6 atoms in bond
   @param atomCoordinates  atom coordinates
   @param parameters       4 parameters: parameters[0] = k
                                         parameters[1] = ideal bond angle in radians
                                         parameters[2] = multiplicity
                                         parameters[3] = cutoff
   @param forces           force array (forces added to current values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCutoffProperDihedralBond::calculateBondIxn(vector<int>& atomIndices,
                                                   vector<Vec3>& atomCoordinates,
                                                   vector<double>& parameters,
                                                   vector<Vec3>& forces,
                                                   double* totalEnergy, double* energyParamDerivs) {
    double deltaR[3][ReferenceForce::LastDeltaRIndex];

    double crossProductMemory[6];

    // ---------------------------------------------------------------------------------------

    // get deltaR, R2, and R between three pairs of atoms: [j,i], [j,k], [l,k]
    int atomAIndex = atomIndices[0];
    int atomBIndex = atomIndices[1];
    int atomCIndex = atomIndices[2];
    int atomDIndex = atomIndices[3];
    if (usePeriodic) {
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], boxVectors, deltaR[0]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], boxVectors, deltaR[1]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], boxVectors, deltaR[2]);
    }
    else {
        ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0]);
        ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], deltaR[1]);
        ReferenceForce::getDeltaR(atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], deltaR[2]);
    }

    double dotDihedral;
    double signOfAngle;
    int hasREntry = 1;

    // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

    double* crossProduct[2];
    crossProduct[0]           = crossProductMemory;
    crossProduct[1]           = crossProductMemory + 3;

    // get dihedral angle

    double dihedralAngle  =  getDihedralAngleBetweenThreeVectors(deltaR[0], deltaR[1], deltaR[2],
                                                                 crossProduct, &dotDihedral, deltaR[0],
                                                                 &signOfAngle, hasREntry);

    // evaluate delta angle, dE/d(angle)

    double deltaAngle     = parameters[2]*dihedralAngle - parameters[1];
    double sinDeltaAngle  = SIN(deltaAngle);
    double dEdAngle       = -parameters[0]*parameters[2]*sinDeltaAngle;
    double energy         =  parameters[0]*(1.0 + cos(deltaAngle));
    dEdAngle*=parameters[3];
    energy*=parameters[3];
    // compute force

    double internalF[4][3];
    double forceFactors[4];
    double normCross1         = DOT3(crossProduct[0], crossProduct[0]);
    double normBC             = deltaR[1][ReferenceForce::RIndex];
    forceFactors[0]    = (-dEdAngle*normBC)/normCross1;

    double normCross2         = DOT3(crossProduct[1], crossProduct[1]);
    forceFactors[3]    = (dEdAngle*normBC)/normCross2;

    forceFactors[1]    = DOT3(deltaR[0], deltaR[1]);
    forceFactors[1]   /= deltaR[1][ReferenceForce::R2Index];

    forceFactors[2]    = DOT3(deltaR[2], deltaR[1]);
    forceFactors[2]   /= deltaR[1][ReferenceForce::R2Index];

    for (int ii = 0; ii < 3; ii++) {

        internalF[0][ii] = forceFactors[0]*crossProduct[0][ii];
        internalF[3][ii] = forceFactors[3]*crossProduct[1][ii];

        double s  = forceFactors[1]*internalF[0][ii] - forceFactors[2]*internalF[3][ii];

        internalF[1][ii] = internalF[0][ii] - s;
        internalF[2][ii] = internalF[3][ii] + s;
    }

    // accumulate forces

    for (int ii = 0; ii < 3; ii++) {
        forces[atomAIndex][ii] += internalF[0][ii];
        forces[atomBIndex][ii] -= internalF[1][ii];
        forces[atomCIndex][ii] -= internalF[2][ii];
        forces[atomDIndex][ii] += internalF[3][ii];
    }

    // accumulate energies

    if (totalEnergy != NULL)
        *totalEnergy += energy;
}
