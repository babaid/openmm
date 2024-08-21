#ifndef OPENMM_CUTOFFANGLEFORCE_H_
#define OPENMM_CUTOFFANGLEFORCE_H_

#include "openmm/Force.h"
#include "openmm/Vec3.h"
#include <map>
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an interaction between groups of three particles that varies Cutoffally with the angle
 * between them.  To use it, create a CutoffAngleForce object then call addAngle() once for each angle.  After
 * an angle has been added, you can modify its force field parameters by calling setAngleParameters().  This will
 * have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT CutoffAngleForce : public Force {
public:
    /**
     * Create a CutoffAngleForce.
     */
    CutoffAngleForce();
    /**
     * Get the number of Cutoff bond angle terms in the potential function
     */
    int getNumAngles() const {
        return angles.size();
    }
    /**
     * Add a bond term to the force field.
     *
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param cutoff        the harmonic force constant for the bond, measured in kJ/mol/nm^2
     * @return the index of the bond that was added
     */
    int addBond(int particle1, int particle2, double cutoff);
    /**
     * Get the force field parameters for a bond term.
     *
     * @param[out] particle1 the index of the first particle connected by the bond
     * @param[out] particle2 the index of the second particle connected by the bond
     * @param[out] cutoff        the harmonic force constant for the bond, measured in kJ/mol/nm^2
     */
    void getBondParameter(int& particle1, int& particle2, double& cutoff) const;
    /**
     * Add an angle term to the force field.
     *
     * @param particle1 the index of the first particle forming the angle
     * @param particle2 the index of the second particle forming the angle
     * @param particle3 the index of the third particle forming the angle
     * @param angle     the equilibrium angle, measured in radians
     * @param k         the Cutoff force constant for the angle, measured in kJ/mol/radian^2
     * @return the index of the angle that was added
     */
    int addAngle(int particle1, int particle2, int particle3, double angle, double k);
    /**
     * Get the force field parameters for an angle term.
     *
     * @param      index     the index of the angle for which to get parameters
     * @param[out] particle1 the index of the first particle forming the angle
     * @param[out] particle2 the index of the second particle forming the angle
     * @param[out] particle3 the index of the third particle forming the angle
     * @param[out] angle     the equilibrium angle, measured in radians
     * @param[out] k         the Cutoff force constant for the angle, measured in kJ/mol/radian^2
     */
    void getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& angle, double& k) const;
    /**
     * Set the force field parameters for an angle term.
     *
     * @param index     the index of the angle for which to set parameters
     * @param particle1 the index of the first particle forming the angle
     * @param particle2 the index of the second particle forming the angle
     * @param particle3 the index of the third particle forming the angle
     * @param angle     the equilibrium angle, measured in radians
     * @param k         the Cutoff force constant for the angle, measured in kJ/mol/radian^2
     */
    void setAngleParameters(int index, int particle1, int particle2, int particle3, double angle, double k);
    /**
     * Update the per-angle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setAngleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-angle parameters.  The set of particles involved
     * in a angle cannot be changed, nor can new angles be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Set whether this force should apply periodic boundary conditions when calculating displacements.
     * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;protected:
    ForceImpl* createImpl() const;
private:
    class BondInfo;
    class AngleInfo;
    std::vector<AngleInfo> angles;
    std::vector<BondInfo> cutoffbond;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about an angle.
 * @private
 */
class CutoffAngleForce::AngleInfo {
public:
    int particle1, particle2, particle3;
    double angle, k;
    AngleInfo() {
        particle1 = particle2 = particle3 = -1;
        angle = k = 0.0;
    }
    AngleInfo(int particle1, int particle2, int particle3, double angle, double k) :
        particle1(particle1), particle2(particle2), particle3(particle3), angle(angle), k(k) {
    }
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class CutoffAngleForce::BondInfo {
public:
    int particle1, particle2;
    double cutoff;
    BondInfo() {
        particle1 = particle2 = -1;
        cutoff = 0.0;
    }
    BondInfo(int particle1, int particle2, double cutoff) :
        particle1(particle1), particle2(particle2),  cutoff(cutoff) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CutoffANGLEFORCE_H_*/
