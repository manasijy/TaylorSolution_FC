phi1 = degtorad(euler_angle(1));
phi = degtorad(euler_angle(2));
phi2 = degtorad(euler_angle(3));

g11 = ((cos(phi1)*cos(phi2))-sin(phi1)*sin(phi2)*cos(phi));
g12 = (sin(phi1)*cos(phi2))+(cos(phi1)*sin(phi2)*cos(phi));
g13 = sin(phi2)*sin(phi);

g21 = (-cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(phi));
g22 = (-sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(phi));
g23 = cos(phi2)*sin(phi);

g31 = sin(phi1)*sin(phi);
g32 = (-cos(phi1)*sin(phi));
g33 = cos(phi);

g = [g11 g12 g13;g21 g22 g23;g31 g32 g33];
