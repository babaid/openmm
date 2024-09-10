real3 sv = make_real3(pos5.x-pos4.x, pos5.y-pos4.y, pos5.z-pos4.z);
real norm_fact_sv = SQRT(sv.x*sv.x + sv.y*sv.y + sv.z*sv.z); 
float3 angleParams = PARAMS[index];
real deltaIdeal = theta-angleParams.x;
real cutoff_factor = 1.0f/(1.0f + __expf(5.0f*(norm_fact_sv - angleParams.z)));
energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal*cutoff_factor;
real dEdAngle = angleParams.y*deltaIdeal*cutoff_factor;