function [DCM313] = DCM313(RAAN_in, inc_in, omega_in)

%Returns DCM for 313 Euler angles (Right ascension of ascending node,
%inclination angle, argument of periapsis)

    %Angles come in as degrees (convention)
    
    RAAN = deg2rad(RAAN_in);
    inc = deg2rad(inc_in);
    omega = deg2rad(omega_in);

    DCM313 = zeros(3,3);
    
    DCM313(1,1) = cos(omega)*cos(RAAN) - sin(omega)*cos(inc)*sin(RAAN);
    DCM313(1,2) = cos(omega)*sin(RAAN) + sin(omega)*cos(inc)*cos(RAAN);
    DCM313(1,3) = sin(omega)*sin(inc);
    
    DCM313(2,1) = -sin(omega)*cos(RAAN) - cos(omega)*cos(inc)*sin(RAAN);
    DCM313(2,2) = -sin(omega)*sin(RAAN) + cos(omega)*cos(inc)*cos(RAAN);
    DCM313(2,3) = cos(omega)*sin(inc);
    
    DCM313(3,1) = sin(inc)*sin(RAAN);
    DCM313(3,2) = -sin(inc)*cos(RAAN);
    DCM313(3,3) = cos(inc);

end