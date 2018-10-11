function [sys] = SRP3BDynamicsParamsSTM(t, X, params)

% State is [R, V, Rs1, Rs2, Rs3, br3, br3dot, muSun, omegaE, Cr]

muEarth = params(1);

% Unpack fixed parameters
PPhi = params(3);
AoverM = params(4);
initialJD = params(5);

% Unpack spacecraft position and station positions
r = X(1:3);
Rs1 = X(7:9);
Rs2 = X(10:12);
Rs3 = X(13:15);

% Unpack estimated parameters from state vector
CR = X(end);
omegaVec = [0; 0; X(end-1)];
muEarth = X(end-2);
muSun = X(end-3);

% Set up dot vector
sys = zeros(length(X),1);

% Spacecraft Velocities
sys(1:3) = X(4:6);

% Get Sun vector
JD = initialJD + (t/86400);

[REarthSun, ~, ~] = Ephem(JD, 3, 'EME2000');
RSun = -REarthSun;

% Spacecraft Accelerations

accelGrav = -muEarth/norm(r)^3*r;
accelSRP = -CR*PPhi*AoverM*(RSun - r)/norm(RSun - r)^3;
accel3BP = muSun*((RSun - r)/norm(RSun - r)^3 - RSun/norm(RSun)^3);

sys(4:6) = accelGrav + accelSRP + accel3BP;

% Station Velocities

sys(7:9) = cross3d(omegaVec, Rs1);
sys(10:12) = cross3d(omegaVec, Rs2);
sys(13:15) = cross3d(omegaVec, Rs3);

% STM
STM = reshape(X(26:end),25,[]);
A = AFullParams(X(1:25),RSun,params);

STMDot = A*STM;

sys(26:end) = STMDot(:);

end