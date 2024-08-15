#include "openmm/serialization/CutoffAngleForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CutoffAngleForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CutoffAngleForceProxy::CutoffAngleForceProxy() : SerializationProxy("CutoffAngleForce") {
}

void CutoffAngleForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const CutoffAngleForce& force = *reinterpret_cast<const CutoffAngleForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
    SerializationNode& cutoffbond = node.createChildNode("CutoffBond");
    int bondingparticle1, bondingparticle2;
    double cutoff;
    force.getBondParameter(bondingparticle1, bondingparticle2, cutoff);
    cutoffbond.createChildNode("Bond").setIntProperty("p1", bondingparticle1).setIntProperty("p2", bondingparticle2).setDoubleProperty("cutoff", cutoff);
    SerializationNode& bonds = node.createChildNode("Angles");
    for (int i = 0; i < force.getNumAngles(); i++) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        bonds.createChildNode("Angle").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setDoubleProperty("a", angle).setDoubleProperty("k", k);
    }
}

void* CutoffAngleForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    CutoffAngleForce* force = new CutoffAngleForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        if (version > 1)
            force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        const SerializationNode& cutoffbond = node.getChildNode("CutoffBond");
        for (auto& bond: cutoffbond.getChildren())
            force->addBond(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getDoubleProperty("cutoff"));
        const SerializationNode& angles = node.getChildNode("Angles");
        for (auto& angle : angles.getChildren())
            force->addAngle(angle.getIntProperty("p1"), angle.getIntProperty("p2"), angle.getIntProperty("p3"), angle.getDoubleProperty("a"), angle.getDoubleProperty("k"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}

