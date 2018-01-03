

function euler = gmat_to_euler(gmat)
phi = real(round(radtodeg(round(acos(gmat(3,3)),4)),4));
phi1 = real(round(radtodeg(round(asin(gmat(3,1)/sin(degtorad(phi))),4)),4));
phi2 = real(round(radtodeg(round(acos(gmat(2,3)/round(sin(degtorad(phi)),3)),4)),4));
euler = round([phi1 phi phi2],3);