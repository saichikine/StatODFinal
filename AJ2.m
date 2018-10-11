function [A] = AJ2(x, y, z, params)

%A = dF/dX, where F = dX/dt, and X = [r v]^T
%so, F = dX/dt = [xdot ydot zdot xddot yddot zddot]^T
%and A = dF/dX = [
%so, A looks like:
%A =  |0    I|
%     |Udd  O|

A = zeros(6,6);

mu = params(1);
J2 = params(2);
R = params(3);

A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

%first row of hessian
A(4,1) = (2*(2*x^2 - y^2 - z^2)*(x^2 + y^2 + z^2)^2 + 3*J2*R^2*(4*x^4 - y^4 + 3*y^2*z^2 + 4*z^4 + 3*x^2*(y^2 - 9*z^2)))*mu/(2*(x^2 + y^2 + z^2)^(9/2));
A(4,2) = (3*mu*x*y*(5*J2*R^2*(x^2+y^2-6*z^2) + 2*(x^2+y^2+z^2)^2))/(2*(x^2+y^2+z^2)^(9/2));
A(4,3) = (3*mu*x*z*(5*J2*R^2*(3*(x^2+y^2) - 4*z^2) + 2*(x^2+y^2+z^2)^2))/(2*(x^2+y^2+z^2)^(9/2));

%second row of hessian
A(5,1) = A(4,2);
A(5,2) = -(2*(x^2 - 2*y^2 + z^2)*(x^2 + y^2 + z^2)^2 + 3*J2*R^2*(x^4 - 4*y^4 + 27*y^2*z^2 - 4*z^4 - 3*x^2*(y^2 + z^2)))*mu/(2*(x^2 + y^2 + z^2)^(9/2));
A(5,3) = (3*mu*y*z*(5*J2*R^2*(3*(x^2+y^2)-4*z^2) + 2*(x^2+y^2+z^2)^2))/(2*(x^2+y^2+z^2)^(9/2));

%third row of hessian
A(6,1) = A(4,3);
A(6,2) = A(5,3);
A(6,3) = -(2*(x^2 + y^2 - 2*z^2)*(x^2 + y^2 + z^2)^2 + 3*J2*R^2*(3*(x^2 + y^2)^2 - 24*(x^2 + y^2)*z^2 + 8*z^4))*mu/(2*(x^2 + y^2 + z^2)^(9/2));

end