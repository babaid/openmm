
#ifndef __ReferenceCutoffProperDihedralBond_H__
#define __ReferenceCutoffProperDihedralBond_H__

#include "ReferenceBondIxn.h"

namespace OpenMM {

    class OPENMM_EXPORT ReferenceCutoffProperDihedralBond : public ReferenceBondIxn {

    private:

        bool usePeriodic;
        Vec3 boxVectors[3];

    public:

        /**---------------------------------------------------------------------------------------

           Constructor

           --------------------------------------------------------------------------------------- */

        ReferenceCutoffProperDihedralBond();

        /**---------------------------------------------------------------------------------------

           Destructor

           --------------------------------------------------------------------------------------- */

        ~ReferenceCutoffProperDihedralBond();

        /**---------------------------------------------------------------------------------------

          Set the force to use periodic boundary conditions.

          @param vectors    the vectors defining the periodic box

          --------------------------------------------------------------------------------------- */

        void setPeriodic(OpenMM::Vec3* vectors);

        /**---------------------------------------------------------------------------------------

          Calculate proper dihedral bond ixn

           @param atomIndices      atom indices of 4 atoms in bond
           @param atomCoordinates  atom coordinates
           @param parameters       3 parameters: parameters[0] = k
                                                 parameters[1] = ideal bond angle in radians
                                                 parameters[2] = multiplicity
           @param forces           force array (forces added to current values)
           @param totalEnergy      if not null, the energy will be added to this

           --------------------------------------------------------------------------------------- */

        void calculateBondIxn(std::vector<int>& atomIndices, std::vector<OpenMM::Vec3>& atomCoordinates,
                              std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces,
                              double* totalEnergy, double* energyParamDerivs);

    };

} // namespace OpenMM

#endif // __ReferenceProperDihedralBond_H__
