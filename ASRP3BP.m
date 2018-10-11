function [A] = ASRP3BP(x, y, z, CR, RSun, params)

muEarth = params(1);
muSun = params(2);
PPhi = params(3);
AoverM = params(4);

A = zeros(7,7);

A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

r = [x; y; z]; % spacecraft position vector

thingX = (norm(RSun-r)^(3)*(-1) - (RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
thingY = (norm(RSun-r)^(3)*(-1) - (RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
thingZ = (norm(RSun-r)^(3)*(-1) - (RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);

A(4,1) = -muEarth*(norm(r)^(3) - x*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingX;
A(4,2) = -muEarth*(-x*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
A(4,3) = -muEarth*(-x*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(1) - x)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);
A(4,7) = -PPhi*AoverM*(RSun(1) - x)/(norm(RSun - r)^(3));

A(5,1) = -muEarth*(-y*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
A(5,2) = -muEarth*(norm(r)^(3) - y*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingY;
A(5,3) = -muEarth*(-y*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(2) - y)*(3/2)*norm(RSun-r)*2*(RSun(3) - z)*(-1))/(norm(RSun-r)^6);
A(5,7) = -PPhi*AoverM*(RSun(2) - y)/(norm(RSun - r)^(3));

A(6,1) = -muEarth*(-z*(3/2)*norm(r)*2*x)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(1) - x)*(-1))/(norm(RSun-r)^6);
A(6,2) = -muEarth*(-z*(3/2)*norm(r)*2*y)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*(-(RSun(3) - z)*(3/2)*norm(RSun-r)*2*(RSun(2) - y)*(-1))/(norm(RSun-r)^6);
A(6,3) = -muEarth*(norm(r)^(3) - z*(3/2)*norm(r)*2*z)/(norm(r)^6) + (muSun - CR*PPhi*AoverM)*thingZ;
A(6,7) = -PPhi*AoverM*(RSun(3) - z)/(norm(RSun - r)^(3));