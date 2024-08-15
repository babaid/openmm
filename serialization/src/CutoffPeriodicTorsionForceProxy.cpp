#include "openmm/serialization/CutoffPeriodicTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CutoffPeriodicTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CutoffPeriodicTorsionForceProxy::CutoffPeriodicTorsionForceProxy() : SerializationProxy("CutoffPeriodicTorsionForce") {
}

void CutoffPeriodicTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const CutoffPeriodicTorsionForce& force = *reinterpret_cast<const CutoffPeriodicTorsionForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());

    SerializationNode& cutoffbond = node.createChildNode("CutoffBond");
    int bondingparticle1, bondingparticle2;
    double cutoff;
    force.getBondParameter(bondingparticle1, bondingparticle2, cutoff);
    cutoffbond.createChildNode("Bond").setIntProperty("p1", bondingparticle1).setIntProperty("p2", bondingparticle2).setDoubleProperty("cutoff", cutoff);

    SerializationNode& torsions = node.createChildNode("Torsions");
    for (int i = 0; i < force.getNumTorsions(); i++) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        torsions.createChildNode("Torsion").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setIntProperty("p4", particle4).setIntProperty("periodicity", periodicity).setDoubleProperty("phase", phase).setDoubleProperty("k", k);
    }
}

void* CutoffPeriodicTorsionForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    CutoffPeriodicTorsionForce* force = new CutoffPeriodicTorsionForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        if (version > 1)
            force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        const SerializationNode& cutoffbond = node.getChildNode("CutoffBond");
        for (auto& bond: cutoffbond.getChildren())
            force->addBond(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getDoubleProperty("cutoff"));
        const SerializationNode& torsions = node.getChildNode("Torsions");
        for (auto& torsion : torsions.getChildren())
            force->addTorsion(torsion.getIntProperty("p1"), torsion.getIntProperty("p2"), torsion.getIntProperty("p3"), torsion.getIntProperty("p4"),
                              torsion.getIntProperty("periodicity"), torsion.getDoubleProperty("phase"), torsion.getDoubleProperty("k"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
