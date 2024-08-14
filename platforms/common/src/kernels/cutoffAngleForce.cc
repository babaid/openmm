real3 v2 = make_real3(pos4.x-pos5.x, pos4.y-pos5.y, pos4.z-pos5.z);
real dist2 = v2.x*v2.x + v2.y*v2.y + v2.z*v2.z;
real dist = SQRT(dist2);
float3 angleParams = PARAMS[index];
real deltaIdeal = theta-angleParams.x;
energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal*(1/(1+exp(dist - angleParams.z)));
real dEdAngle = angleParams.y*deltaIdeal*(1/(1+exp(dist - angleParams.z)));