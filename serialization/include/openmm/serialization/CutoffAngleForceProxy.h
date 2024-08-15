#ifndef OPENMM_CUTOFFANGLEFORCE_PROXY_H_
#define OPENMM_CUTOFFANGLEFORCE_PROXY_H_


#include "openmm/internal/windowsExport.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

/**
 * This is a proxy for serializing CutoffAngleForce objects.
 */

    class OPENMM_EXPORT CutoffAngleForceProxy : public SerializationProxy {
    public:
        CutoffAngleForceProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

} // namespace OpenMM

#endif /*OPENMM_HARMONICANGLEFORCE_PROXY_H_*/
