#ifndef OPENMM_CUTOFFPERIODICTORSIONFORCEIMPL_H_
#define OPENMM_CUTOFFPERIODICTORSIONFORCEIMPL_H_


#include "openmm/internal/ForceImpl.h"
#include "../CutoffPeriodicTorsionForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of PeriodicTorsionForce.
 */

    class CutoffPeriodicTorsionForceImpl : public ForceImpl {
    public:
        CutoffPeriodicTorsionForceImpl(const CutoffPeriodicTorsionForce& owner);
        ~CutoffPeriodicTorsionForceImpl();
        void initialize(ContextImpl& context);
        const CutoffPeriodicTorsionForce& getOwner() const {
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
        const CutoffPeriodicTorsionForce& owner;
        Kernel kernel;
    };

} // namespace OpenMM

#endif /*OPENMM_CUTOFFPERIODICTORSIONFORCEIMPL_H_*/
