float3 angleParams = PARAMS[index];
real deltaIdeal = theta-angleParams.x;
energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal*angleParams.z;
real dEdAngle = angleParams.y*deltaIdeal*angleParams.z;