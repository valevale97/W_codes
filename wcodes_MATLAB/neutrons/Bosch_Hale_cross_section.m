function sigma = Bosch_Hale_cross_section(E)
% This script calculates the D-D fusion cross section using the
% approximation published by Bosch and Hale in NF 1992. E is the energy in the
% center-of-mass frame in units of keV. For a stationary target ion, the
% center-of-mass energy is simply half the energy of the fast ion in the
% lab frame. The cross section is given in units of mb (milli-barns). 
% 1 mb = 10^-31 m^2

B_G = 31.3970;

A1 = 5.3701e4;
A2 = 3.3027e2;
A3 = -1.2706e-1;
A4 = 2.9327e-5;
A5 = -2.5151e-9;

S = A1 + E.*(A2 + E.*(A3 + E.*(A4 + E.*A5)));

sigma = S./(E.*exp(B_G./sqrt(E)));

end