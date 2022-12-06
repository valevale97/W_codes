function sigmadiff = diffCrossSection(E1, phi)
% Compute the DD differential cross section in the lab frame, where one of
% the deuterons is assumed to be at rest. 'E1' is the energy of the
% incomming particle in MeV and 'phi' is the emission angle in degrees, 
% relative to the incident velocity vector. The cross section is returned
% in units of millibarn/sr.

% Constants
c   = 299792458;              % speed of light [m/s]
MeV = 1.60217e-13;            % [J/MeV]
m1  = 3.34358e-27;            % deuteron mass [kg]
m2  = m1;
m3  = 1.67493e-27;            % neutron mass [kg]
m4  = 5.00641e-27;            % He-3 mass [kg]
Q   = 3.26893;                % fusion Q-value [MeV]

% LAB quanteties
v1   = sqrt(2*E1*MeV/m1);     % speed of incomming particle [m/s]

% CMS quanteties
K      = E1 /2.0;                % total kinetic energy [MeV]
vcm    = v1 / 2.0;               % speed of CMS [m/s]
E3cm   = m4/(m3+m4) * (K + Q);   % energy of the neutron [MeV]
v3cm   = sqrt(2*E3cm*MeV/m3);    % speed of the neutron [m/s]

x = vcm ./ v3cm;
cosTheta = cosd(phi) .* sqrt(1 - x.^2.*sind(phi).^2) - x.*sind(phi).^2;     % emission angle in CMS

cosTheta(cosTheta<-1) = -1;       % correct for numerical effects
cosTheta(cosTheta>1) = 1;

% Angular dependence from Lagendre parameterization
load legendrecoeff.mat
A = interp1(legendrecoeff(:,1), legendrecoeff(:,2:end), K, 'linear', 0);
P = zeros(size(A));

for i = 1:6
    temp = legendre(2*i-2,cosTheta)';
    P(:,i) = temp(:,1);
end

C = sum(A .* P, 2) ./ (4*pi * A(:,1));

% Bosch&Hale parameterization of the total cross section
B_G = 31.3970;

A1  = 5.3701e4;
A2  = 3.3027e2;
A3  = -1.2706e-1;
A4  = 2.9327e-5;
A5  = -2.5151e-9;

E = K * 1000;   % convert to keV
S   = A1+E.*(A2+E.*(A3+E.*(A4+E.*A5)));
sigma = S./(E.*exp(B_G./sqrt(E)));

% Multiply the angular part with the total cross section
sigmadiff = C' .* sigma;    % diff cross section in CMS

% Convert to the LAB frame, assuming the target particle to be at rest.
J = (1 + x.^2 + 2*x.*cosTheta).^1.5 ./ (1 + x.*cosTheta);                       % Jacobian for the transformation
sigmadiff = sigmadiff .* J;

end