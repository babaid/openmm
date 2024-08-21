#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "../../plugins/semibondforces/openmmapi/include/openmm/CutoffAngleForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "../../plugins/semibondforces/openmmapi/include/openmm/internal/CutoffAngleForceImpl.h"

using namespace OpenMM;

CutoffAngleForce::CutoffAngleForce() : usePeriodic(false) {
}
int CutoffAngleForce::addBond(int particle1, int particle2, double cutoff) {
    cutoffbond.push_back(BondInfo(particle1, particle2, cutoff));
    return 0;
}

void CutoffAngleForce::getBondParameter(int &particle1, int &particle2, double &cutoff) const {
    particle1 = cutoffbond[0].particle1;
    particle2 = cutoffbond[0].particle2;
    cutoff = cutoffbond[0].cutoff;
}
int CutoffAngleForce::addAngle(int particle1, int particle2, int particle3, double angle, double k) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, angle, k));
    return angles.size()-1;
}

void CutoffAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& angle, double& k) const {
    ASSERT_VALID_INDEX(index, angles);
    particle1 = angles[index].particle1;
    particle2 = angles[index].particle2;
    particle3 = angles[index].particle3;
    angle = angles[index].angle;
    k = angles[index].k;
}

void CutoffAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, double angle, double k) {
    ASSERT_VALID_INDEX(index, angles);
    angles[index].particle1 = particle1;
    angles[index].particle2 = particle2;
    angles[index].particle3 = particle3;
    angles[index].angle = angle;
    angles[index].k = k;
}

ForceImpl* CutoffAngleForce::createImpl() const {
    return new CutoffAngleForceImpl(*this);
}

void CutoffAngleForce::updateParametersInContext(Context& context) {
    dynamic_cast<CutoffAngleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CutoffAngleForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CutoffAngleForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
