#ifndef OPENMM_CUTOFFPERIODICTORSIONFORCE_PROXY_H_
#define OPENMM_CUTOFFPERIODICTORSIONFORCE_PROXY_H_

#include "openmm/internal/windowsExport.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

/**
 * This is a proxy for serializing PeriodicTorsionForce objects.
 */

    class OPENMM_EXPORT CutoffPeriodicTorsionForceProxy : public SerializationProxy {
    public:
    CutoffPeriodicTorsionForceProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMM

#endif /*OPENMM_PERIODICTORSIONFORCE_PROXY_H_*/
