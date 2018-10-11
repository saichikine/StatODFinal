function [sys] = SRP3BDynamics_s(t, X, RSun, params)

% State is [R, V, Rs, Cr]

muEarth = params(1);
muSun = params(2);
PPhi = params(3);
AoverM = params(4);
initialJD = params(5);

sys = zeros(length(X),1);

% Spacecraft Velocities
sys(1:3) = X(4:6);

x = X(1);
y = X(2);
z = X(3);

Rs1 = X(7:9);
Rs2 = X(10:12);
Rs3 = X(13:15);

r = [x; y; z];

% Get Sun vector
% JD = initialJD + (t/86400);
% 
% [REarthSun, ~, ~] = Ephem(JD, 3, 'EME2000');
% RSun = -REarthSun;

% Spacecraft Accelerations

CR = X(end);

accelGrav = -muEarth/norm(r)^3*r;
accelSRP = -CR*PPhi*AoverM*(RSun - r)/norm(RSun - r)^3;
accel3BP = muSun*((RSun - r)/norm(RSun - r)^3 - RSun/norm(RSun)^3);

sys(4:6) = accelGrav + accelSRP + accel3BP;

% Station Velocities

omegaVec = [0; 0; 7.29211585275553e-5];
sys(7:9) = cross(omegaVec, Rs1);
sys(10:12) = cross(omegaVec, Rs2);
sys(13:15) = cross(omegaVec, Rs3);

end