function [sys] = J3DynamicsSTM(t, X, params)

mu = params(1);
J2 = params(2);
R = params(3);
J3 = params(4);

sys = zeros(42,1);

%Velocities
sys(1:3) = X(4:6);

x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);

sys(4:6) = [(-1/2).*x.*(x.^2+y.^2+z.^2).^(-9/2).*(2.*(x.^2+y.^2).^3+15.*J3.* ...
  R.^3.*(x.^2+y.^2).*z+6.*(x.^2+y.^2).^2.*z.^2+(-20).*J3.*R.^3.* ...
 z.^3+6.*(x.^2+y.^2).*z.^4+2.*z.^6+3.*J2.*R.^2.*(x.^2+y.^2+(-4).* ...
  z.^2).*(x.^2+y.^2+z.^2)).*mu,(-1/2).*y.*(x.^2+y.^2+z.^2).^(-9/2).*( ...
  2.*(x.^2+y.^2).^3+15.*J3.*R.^3.*(x.^2+y.^2).*z+6.*(x.^2+y.^2).^2.* ...
z.^2+(-20).*J3.*R.^3.*z.^3+6.*(x.^2+y.^2).*z.^4+2.*z.^6+3.*J2.* ...
R.^2.*(x.^2+y.^2+(-4).*z.^2).*(x.^2+y.^2+z.^2)).*mu,(-1/2).*(x.^2+ ...
  y.^2+z.^2).^(-9/2).*(J3.*R.^3.*((-3).*(x.^2+y.^2).^2+24.*(x.^2+ ...
  y.^2).*z.^2+(-8).*z.^4)+z.*(x.^2+y.^2+z.^2).*(3.*J2.*R.^2.*(3.*( ...
  x.^2+y.^2)+(-2).*z.^2)+2.*(x.^2+y.^2+z.^2).^2)).*mu]';

%Accelerations
% sys(4) = -mu*x*(3*J2*R^2*(x^2+y^2-4*z^2) + 2*(x^2+y^2+z^2)^2)/(2*(x^2+y^2+z^2)^(7/2));
% sys(5) = -mu*y*(3*J2*R^2*(x^2+y^2-4*z^2) + 2*(x^2+y^2+z^2)^2)/(2*(x^2+y^2+z^2)^(7/2));
% sys(6) = -mu*z*(3*J2*R^2*(3*(x^2+y^2) - 2*z^2) + 2*(x^2+y^2+z^2)^2)/(2*(x^2+y^2+z^2)^(7/2));

%Get A(t)
At = AJ3(x, y, z, params);

%Get STM from state vector
Phi = reshape(X(7:end), 6, []);

%phidot = A*phi
Phi_dot = At*Phi;

%reshape phidot to a column vector
Phi_dot_column = reshape(Phi_dot, [], 1);

%Append phidot column vector to state vector
sys(7:end) = Phi_dot_column;
end