function [A] = AFullParams(X, RSun, params)

PPhi = params(3);
AoverM = params(4);

% Unpack estimated parameters from state vector
CR = X(end);
omegaVec = [0; 0; X(end-1)];
muEarth = X(end-2);
muSun = X(end-3);

L = length(X);

A = zeros(L,L);

A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

r = X(1:3); % spacecraft position vector
omegaTilde = [0 -omegaVec(3) 0; omegaVec(3) 0 0; 0 0 0];

x = r(1);
y = r(2);
z = r(3);

thingX = (norm(RSun-r)^(3)*(-1) - (RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
thingY = (norm(RSun-r)^(3)*(-1) - (RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
thingZ = (norm(RSun-r)^(3)*(-1) - (RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);

A(4,1) = -muEarth*(norm(r)^(3) - x*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingX;
A(4,2) = -muEarth*(-x*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
A(4,3) = -muEarth*(-x*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);
A(4,23) = -1/norm(r)^3*r(1); %muEarth deriv
A(4,25) = -PPhi*AoverM*(RSun(1) - x)/(norm(RSun - r)^(3)); %CR deriv

A(5,1) = -muEarth*(-y*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
A(5,2) = -muEarth*(norm(r)^(3) - y*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingY;
A(5,3) = -muEarth*(-y*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);
A(5,23) = -1/norm(r)^3*r(2); %muEarth deriv
A(5,25) = -PPhi*AoverM*(RSun(2) - y)/(norm(RSun - r)^(3)); %CR deriv

A(6,1) = -muEarth*(-z*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
A(6,2) = -muEarth*(-z*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
A(6,3) = -muEarth*(norm(r)^(3) - z*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingZ;
A(6,23) = -1/norm(r)^3*r(3); %muEarth deriv
A(6,25) = -PPhi*AoverM*(RSun(3) - z)/(norm(RSun - r)^(3)); %CR Deriv

%muSun derivs
A(4:6,22) = ((RSun - r)/norm(RSun - r)^3 - RSun/norm(RSun)^3);

%derivs of station velocites wrt to omega
A(7:15,7:15) = blkdiag(omegaTilde, omegaTilde, omegaTilde);


