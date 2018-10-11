function [sys] = SRP3BDynamicsParamsDMC(t, X, RSun, params)

% State is [R, V, Rs1, Rs2, Rs3, br3, br3dot, muSun, omegaE, Cr]

% Unpack fixed parameters
PPhi = params(3);
AoverM = params(4);
tau = params(6);

% Unpack spacecraft position and station positions
r = X(1:3);
Rs1 = X(10:12);
Rs2 = X(13:15);
Rs3 = X(16:18);

% Unpack estimated parameters from state vector
CR = X(25+3);
omegaVec = [0; 0; X(24+3)];
muEarth = X(23+3);
muSun = X(22+3);

% Set up dot vector
sys = zeros(length(X),1);

% Spacecraft Velocities
sys(1:3) = X(4:6);

% Spacecraft Accelerations

accelDMC = X(7:9);

accelGrav = -muEarth/norm(r)^3*r;
accelSRP = -CR*PPhi*AoverM*(RSun - r)/norm(RSun - r)^3;
accel3BP = muSun*((RSun - r)/norm(RSun - r)^3 - RSun/norm(RSun)^3);

% threshold
%accelDMC = (1e-17*(abs(accelDMC) < 1e-17)).*sign(accelDMC) + (accelDMC .* (abs(accelDMC) >= 1e-17));

% DMC Acceleration Rates
sys(7:9) = -1/tau*accelDMC;

sys(4:6) = accelGrav + accelSRP + accel3BP + accelDMC;

% Station Velocities

sys(10:12) = cross3d(omegaVec, Rs1);
sys(13:15) = cross3d(omegaVec, Rs2);
sys(16:18) = cross3d(omegaVec, Rs3);

end