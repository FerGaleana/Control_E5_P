function REF = DOF_PLANAR_P2_REFu( t )
N = 5;
ref1 = (pi/3)*sin(N*t);
ref2 = (pi/3)*cos(N*t);

REF = [ref1 , ref2];

end