real3 vdist = make_real3(pos5.x-pos6.x, pos5.y-pos6.y, pos5.z-pos6.z);
real dist2 = vdist.x*vdist.x + vdist.y*vdist.y + vdist.z*vdist.z;
real dist = SQRT(dist2);

float4 torsionParams = PARAMS[index];
real deltaAngle = torsionParams.z*theta-torsionParams.y;
energy += torsionParams.x*(1.0f+COS(deltaAngle))*torsionParams.w;
real sinDeltaAngle = SIN(deltaAngle);
real dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle*torsionParams.w;
