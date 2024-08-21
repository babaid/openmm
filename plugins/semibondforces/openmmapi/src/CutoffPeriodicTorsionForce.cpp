#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "../../plugins/semibondforces/openmmapi/include/openmm/CutoffPeriodicTorsionForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "../../plugins/semibondforces/openmmapi/include/openmm/internal/CutoffPeriodicTorsionForceImpl.h"

using namespace OpenMM;

CutoffPeriodicTorsionForce::CutoffPeriodicTorsionForce() : usePeriodic(false) {
}
int CutoffPeriodicTorsionForce::addBond(int particle1, int particle2, double cutoff) {
    cutoffbond.push_back(BondInfo(particle1, particle2, cutoff));
    return 0;
}

void CutoffPeriodicTorsionForce::getBondParameter(int &particle1, int &particle2, double &cutoff) const {
    particle1 = cutoffbond[0].particle1;
    particle2 = cutoffbond[0].particle2;
    cutoff = cutoffbond[0].cutoff;
}

int CutoffPeriodicTorsionForce::addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    periodicTorsions.push_back(PeriodicTorsionInfo(particle1, particle2, particle3, particle4, periodicity, phase, k));
    return periodicTorsions.size()-1;
}

void CutoffPeriodicTorsionForce::getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& periodicity, double& phase, double& k) const {
    ASSERT_VALID_INDEX(index, periodicTorsions);
    particle1 = periodicTorsions[index].particle1;
    particle2 = periodicTorsions[index].particle2;
    particle3 = periodicTorsions[index].particle3;
    particle4 = periodicTorsions[index].particle4;
    periodicity = periodicTorsions[index].periodicity;
    phase = periodicTorsions[index].phase;
    k = periodicTorsions[index].k;
}

void CutoffPeriodicTorsionForce::setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    ASSERT_VALID_INDEX(index, periodicTorsions);
    periodicTorsions[index].particle1 = particle1;
    periodicTorsions[index].particle2 = particle2;
    periodicTorsions[index].particle3 = particle3;
    periodicTorsions[index].particle4 = particle4;
    periodicTorsions[index].periodicity = periodicity;
    periodicTorsions[index].phase = phase;
    periodicTorsions[index].k = k;
}

ForceImpl* CutoffPeriodicTorsionForce::createImpl() const {
    return new CutoffPeriodicTorsionForceImpl(*this);
}

void CutoffPeriodicTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<CutoffPeriodicTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CutoffPeriodicTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CutoffPeriodicTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
