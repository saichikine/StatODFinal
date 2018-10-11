function [R V] = COEstoRV(a, e, i, RAAN, omega, nu, mu)
    p = a*(1-e^2); %semi-latus rectum, [km]
    r = p/(1+e*cos(nu)); %orbit radius, [km]

    h = sqrt(mu*a*(1-e^2)); %angular momentum

    x = r*(cos(RAAN)*cos(omega+nu) - sin(RAAN)*sin(omega+nu)*cos(i)); %x-position, [km]
    y = r*(sin(RAAN)*cos(omega+nu) + cos(RAAN)*sin(omega+nu)*cos(i)); %y-position, [km]
    z = r*(sin(i)*sin(omega+nu)); %z-position, [km]

    xdot = x*h*e/(r*p)*sin(nu) - h/r*(cos(RAAN)*sin(omega+nu) + sin(RAAN)*cos(omega+nu)*cos(i)); %x-velocity, [km/s]
    ydot = y*h*e/(r*p)*sin(nu) - h/r*(sin(RAAN)*sin(omega+nu) - cos(RAAN)*cos(omega+nu)*cos(i)); %y-velocity, [km/s]
    zdot = z*h*e/(r*p)*sin(nu) + h/r*sin(i)*cos(omega+nu);
    
    R = [x; y; z];
    V = [xdot; ydot; zdot];
    
    %[R V] = [[x y z] [xdot ydot zdot]]'; %cartesian state vector
end