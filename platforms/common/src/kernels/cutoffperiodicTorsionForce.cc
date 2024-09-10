real3 sv = make_real3(pos6.x-pos5.x, pos6.y-pos5.y, pos6.z-pos5.z);
real norm_fact_sv = SQRT(sv.x*sv.x + sv.y*sv.y + sv.z*sv.z); 
float4 torsionParams = PARAMS[index];
real deltaAngle = torsionParams.z*theta-torsionParams.y;
real cutoff_factor = 1.0f/(1.0f + __expf(5.0f*(norm_fact_sv -  torsionParams.w)));
energy += torsionParams.x*(1.0f+COS(deltaAngle))*cutoff_factor;
real sinDeltaAngle = SIN(deltaAngle);
real dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle*cutoff_factor;
