#ifndef OPENMM_CUTOFFPERIODICTORSIONFORCE_H_
#define OPENMM_CUTOFFPERIODICTORSIONFORCE_H_


#include "Force.h"
#include "Vec3.h"
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an interaction between groups of four particles that varies periodically with the torsion angle
 * between them.  To use it, create a PeriodicTorsionForce object then call addTorsion() once for each torsion.  After
 * a torsion has been added, you can modify its force field parameters by calling setTorsionParameters().  This will
 * have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT CutoffPeriodicTorsionForce : public Force {
public:
    /**
     * Create a PeriodicTorsionForce.
     */
    CutoffPeriodicTorsionForce();
    /**
     * Get the number of periodic torsion terms in the potential function
     */
    int getNumTorsions() const {
        return periodicTorsions.size();
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
     * Add a periodic torsion term to the force field.
     *
     * @param particle1    the index of the first particle forming the torsion
     * @param particle2    the index of the second particle forming the torsion
     * @param particle3    the index of the third particle forming the torsion
     * @param particle4    the index of the fourth particle forming the torsion
     * @param periodicity  the periodicity of the torsion
     * @param phase        the phase offset of the torsion, measured in radians
     * @param k            the force constant for the torsion
     * @return the index of the torsion that was added
     */
    int addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k);
    /**
     * Get the force field parameters for a periodic torsion term.
     *
     * @param index             the index of the torsion for which to get parameters
     * @param[out] particle1    the index of the first particle forming the torsion
     * @param[out] particle2    the index of the second particle forming the torsion
     * @param[out] particle3    the index of the third particle forming the torsion
     * @param[out] particle4    the index of the fourth particle forming the torsion
     * @param[out] periodicity  the periodicity of the torsion
     * @param[out] phase        the phase offset of the torsion, measured in radians
     * @param[out] k            the force constant for the torsion
     */
    void getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& periodicity, double& phase, double& k) const;
    /**
     * Set the force field parameters for a periodic torsion term.
     *
     * @param index        the index of the torsion for which to set parameters
     * @param particle1    the index of the first particle forming the torsion
     * @param particle2    the index of the second particle forming the torsion
     * @param particle3    the index of the third particle forming the torsion
     * @param particle4    the index of the fourth particle forming the torsion
     * @param periodicity  the periodicity of the torsion
     * @param phase        the phase offset of the torsion, measured in radians
     * @param k            the force constant for the torsion
     */
    void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k);
    /**
     * Update the per-torsion parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setTorsionParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-torsion parameters.  The set of particles involved
     * in a torsion cannot be changed, nor can new torsions be added.
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
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class BondInfo;
    class PeriodicTorsionInfo;
    std::vector<PeriodicTorsionInfo> periodicTorsions;
    std::vector<BondInfo> cutoffbond;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a torsion.
 * @private
 */
    class CutoffPeriodicTorsionForce::PeriodicTorsionInfo {
    public:
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        PeriodicTorsionInfo() {
            particle1 = particle2 = particle3 = particle4 = -1;
            periodicity = 1;
            phase = k = 0.0;
        }
        PeriodicTorsionInfo(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) :
                particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), periodicity(periodicity), phase(phase), k(k) {
        }
    };
    /**
 * This is an internal class used to record information about a bond.
 * @private
 */
    class CutoffPeriodicTorsionForce::BondInfo {
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

#endif /*OPENMM_PERIODICTORSIONFORCE_H_*/
