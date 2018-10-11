function [sys] = SRP3BDynamicsParamsSunIn(t, X, params)

% State is [R, V, Rs1, Rs2, Rs3, br3, br3dot, muSun, omegaE, Cr]

% Unpack fixed parameters
muEarth = params(1);
PPhi = params(3);
AoverM = params(4);
initialJD = params(5);

% Unpack spacecraft position and station positions
r = X(1:3);
Rs1 = X(7:9);
Rs2 = X(10:12);
Rs3 = X(13:15);

% Get Sun vector
JD = initialJD + (t/86400);

[REarthSun, ~, ~] = Ephem(JD, 3, 'EME2000');
RSun = -REarthSun;

% Unpack estimated parameters from state vector
CR = X(end);
omegaVec = [0; 0; X(end-1)];
muSun = X(end-3);
muEarth = X(end-2);

% Set up dot vector
sys = zeros(length(X),1);

% Spacecraft Velocities
sys(1:3) = X(4:6);

% Spacecraft Accelerations

accelGrav = -muEarth/norm(r)^3*r;
accelSRP = -CR*PPhi*AoverM*(RSun - r)/norm(RSun - r)^3;
accel3BP = muSun*((RSun - r)/norm(RSun - r)^3 - RSun/norm(RSun)^3);

sys(4:6) = accelGrav + accelSRP + accel3BP;

% Station Velocities

sys(7:9) = cross3d(omegaVec, Rs1);
sys(10:12) = cross3d(omegaVec, Rs2);
sys(13:15) = cross3d(omegaVec, Rs3);

end