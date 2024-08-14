#ifndef __ReferenceCutoffAngleBondIxn_H__
#define __ReferenceCutoffAngleBondIxn_H__

#include "ReferenceBondIxn.h"
#include "openmm/Vec3.h"

namespace OpenMM {

    class OPENMM_EXPORT ReferenceCutoffAngleBondIxn : public ReferenceBondIxn {

    private:

        bool usePeriodic;
        Vec3 boxVectors[3];

    public:

        /**---------------------------------------------------------------------------------------

           Constructor

           --------------------------------------------------------------------------------------- */

        ReferenceCutoffAngleBondIxn();

        /**---------------------------------------------------------------------------------------

           Destructor

           --------------------------------------------------------------------------------------- */

        ~ReferenceCutoffAngleBondIxn();

        /**---------------------------------------------------------------------------------------

          Set the force to use periodic boundary conditions.

          @param vectors    the vectors defining the periodic box

          --------------------------------------------------------------------------------------- */

        void setPeriodic(OpenMM::Vec3* vectors);

        /**---------------------------------------------------------------------------------------

           Get dEdR and energy term for angle bond

           @param  cosine               cosine of angle
           @param  angleParameters      angleParameters: angleParameters[0] = angle in radians
                                                         angleParameters[1] = k (force constant)
           @param  dEdR                 output dEdR
           @param  energyTerm           output energyTerm

           --------------------------------------------------------------------------------------- */

        void getPrefactorsGivenAngleCosine(double cosine, std::vector<double>& angleParameters,
                                           double* dEdR, double* energyTerm) const;

        /**---------------------------------------------------------------------------------------

           Calculate Angle Bond ixn

           @param atomIndices      two bond indices
           @param atomCoordinates  atom coordinates
           @param parameters       parameters: parameters[0] = ideal bond length
                                               parameters[1] = bond k (includes factor of 2)
           @param forces           force array (forces added)
           @param totalEnergy      if not null, the energy will be added to this

           --------------------------------------------------------------------------------------- */

        void calculateBondIxn(std::vector<int>& atomIndices, std::vector<OpenMM::Vec3>& atomCoordinates,
                              std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces,
                              double* totalEnergy, double* energyParamDerivs);


    };

} // namespace OpenMM

#endif // __ReferenceAngleBondIxn_H__
