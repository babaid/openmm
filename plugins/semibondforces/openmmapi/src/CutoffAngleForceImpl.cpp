#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "../../plugins/semibondforces/openmmapi/include/openmm/internal/CutoffAngleForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using namespace std;

CutoffAngleForceImpl::CutoffAngleForceImpl(const CutoffAngleForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

CutoffAngleForceImpl::~CutoffAngleForceImpl() {
}

void CutoffAngleForceImpl::initialize(ContextImpl& context) {
    int bondparticle1, bondaparticle2;
    double cutoff;
    owner.getBondParameter(bondparticle1, bondaparticle2, cutoff);
    const System& system = context.getSystem();
    for (int i = 0; i < owner.getNumAngles(); i++) {
        int particle[3];
        double angle, k;
        owner.getAngleParameters(i, particle[0], particle[1], particle[2], angle, k);
        for (int j = 0; j < 3; j++) {
            if (particle[j] < 0 || particle[j] >= system.getNumParticles()) {
                stringstream msg;
                msg << "CutoffAngleForce: Illegal particle index for an angle: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (angle < 0 || angle > M_PI*1.000001)
            throw OpenMMException("CutoffAngleForce: angle must be between 0 and pi");
    }
    kernel = context.getPlatform().createKernel(CalcCutoffAngleForceKernel::Name(), context);
    kernel.getAs<CalcCutoffAngleForceKernel>().initialize(context.getSystem(), owner);
}

double CutoffAngleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcCutoffAngleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> CutoffAngleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcCutoffAngleForceKernel::Name());
    return names;
}

void CutoffAngleForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCutoffAngleForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
