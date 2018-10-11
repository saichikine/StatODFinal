function BPlaneStuff = BPlaneCompute(state, mu)

    rVec = state(1:3);
    vVec = state(4:6);
    r = norm(rVec);
    v = norm(vVec);
    
    eVec = 1/mu*((v^2 - mu/r)*rVec - dot(rVec,vVec)*vVec);
    e = norm(eVec);
    
    PHat = eVec/e;
    
    hVec = cross(rVec, vVec);
    h = norm(hVec);
    
    WHat = hVec/h;
    QHat = cross(WHat, PHat);
    
    p = h^2/mu;
    
    %a = p/(1-e);
    a = 1/(2/r - v^2/mu);
    
    b = abs(a)*sqrt(e^2 - 1);
    
    SHat = vVec/v;
    NHat = [0; 0; 1];
    THat = cross(SHat, NHat)/norm(cross(SHat, NHat));
    RHat = cross(SHat, THat);
    
    s = -dot(rVec, vVec)/v^2;
    
    B = b*cross(SHat, WHat);
    %B = rVec + s*vVec;
    
    % Find rotation matrix from ECI frame to B-Plane frame (S, T, R) 
    BPlaneDCM = [SHat, THat, RHat]';
    
    nu = acos(dot(rVec/r, PHat));
    f = acosh(1 + v^2/mu*a*(1-e^2)/(1+e*cos(nu)));
    LTOF = mu/v^3*((sinh(f) - f));
    
    BPlaneStuff = {SHat, THat, RHat, BPlaneDCM, B, LTOF};

end