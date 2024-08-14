#ifndef OPENMM_CUTOFANGLEFORCEIMPL_H_
#define OPENMM_CUTOFANGLEFORCEIMPL_H_



#include "ForceImpl.h"
#include "openmm/CutoffAngleForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of CutoffAngleForce.
 */

class CutoffAngleForceImpl : public ForceImpl {
public:
    CutoffAngleForceImpl(const CutoffAngleForce& owner);
    ~CutoffAngleForceImpl();
    void initialize(ContextImpl& context);
    const CutoffAngleForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
private:
    const CutoffAngleForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_CutoffANGLEFORCEIMPL_H_*/
