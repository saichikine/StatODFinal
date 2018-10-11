function BPlaneProblem = EstimateBPlaneParameters(t0, X0, P, dist, params)
    
    RSOI3EventOpts = odeset('AbsTol', 1e-20, 'RelTol', 1e-13, 'Events', @(t,X) EventRSOI3(t,X,dist));

    ICs = [X0; reshape(eye(7), [], 1)];

    [tRSOI3, XRSOI3] = ode45(@(t,X) SRP3BDynamicsSTM(t, X, params), [t0, inf], ICs, RSOI3EventOpts);
    
    % B Plane Stuff
    muEarth = params(1);
    BPlaneStuff = BPlaneCompute(XRSOI3(end,1:6)', muEarth);
    BPlaneDCM = BPlaneStuff{4};
    SHat = BPlaneStuff{1};

    % Integrate to B-Plane crossing

    BPlaneCrossingEventOpts = odeset('AbsTol', 1e-20, 'RelTol', 1e-13, 'Events', @(t,X) EventBPlaneCrossing(t,X,SHat));
    ICs = XRSOI3(end,:);

    [tBPlane, XBPlane] = ode45(@(t,X) SRP3BDynamicsSTM(t, X, params), [tRSOI3(end), inf], ICs, BPlaneCrossingEventOpts);

    STM = reshape(XBPlane(end,8:end),7,7);
    P = STM(1:6,1:6)*P*STM(1:6,1:6)';
    
    PCrossingBPlaneFrame = BPlaneDCM*P(1:3,1:3)*BPlaneDCM';
    PCrossingBPlaneFrameTR = PCrossingBPlaneFrame(2:3,2:3);
    
    BdotR = dot(BPlaneStuff{5}, BPlaneStuff{3});
    BdotT = dot(BPlaneStuff{5}, BPlaneStuff{2});

    BPlaneParameters = [BdotR, BdotT];
    
    BPlaneProblem = {BPlaneParameters, PCrossingBPlaneFrameTR};

end