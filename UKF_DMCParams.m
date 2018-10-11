function UKFStruct = UKF_DMCParams(initStuff, measurements, noiseStuff, dynamicsFunction, measurementFunction, ...
    UKFParams, envParams, varargin)
    
    % Create anonymous function for dynamics
    dynamics = @(t,X,RSun)dynamicsFunction(t,X,RSun,envParams);
    
    function [augSys] = augmentedDynamics(t, X) 
        % This function lets us integrate all the sigma points without
        % having to create a modified dynamics function
        
        augSys = NaN(L*(2*L+1),1);
        
        % Get Sun vector
        JD = envParams(5) + (t/86400);

        [REarthSun, ~, ~] = Ephem(JD, 3, 'EME2000');
        RSun = -REarthSun;
        
        for ind = 1:(2*L+1)
            state = X((L*ind-(L-1)):(L*ind));
            augSys((L*ind-(L-1)):(L*ind)) = dynamics(t, state, RSun);
        end
    end
    
    %% Handle optional arguments
    numVarArgs = length(varargin);
    if numVarArgs > 3
        error('UKF:TooManyInputs', 'requires at most 2 optional arguments');
    end
    
    optArgs = {odeset('RelTol', 1e-13, 'AbsTol', 1e-18), 0};
    optArgs(1:numVarArgs) = varargin;
    [odeOpts, XTruth] = optArgs{:};
    
    if XTruth == 0
        truthFlag = 0;
    else
        truthFlag = 1;
    end
    
    %% Input Conditioning
    X0Ref = initStuff{1};
    P0 = initStuff{2};
    R = initStuff{3};
    startIndex = initStuff{4};
    stopIndex = initStuff{5};
    tau = envParams(6);
    %Q = noiseStuff{1};
    
    alpha = UKFParams(1);
    Beta = UKFParams(2);
    beta = 1/tau;
    
    L = length(X0Ref); % length of state vector
    yLen = 2; % length of measurement vector - don't hard code this
    %measurements = y(startIndex:stopIndex,:);
    N = length(measurements);
    tVec = measurements(:,1);
    %stationStates = stationStatesAll(:,:,startIndex:stopIndex);

    %% Create Result Variables
    
    % States and Covar
    XHist = NaN(L,N);
    XHist(:,1) = X0Ref;
    PHist = NaN(L,L,N);
    PHist(:,:,1) = P0;
    sigmaHist = NaN(L,N);
    sigmaHist(:,1) = sqrt(diag(P0));
    
    % Residuals
    prefitResiduals = NaN(2,N,3);
    postfitResiduals = NaN(2,N,3);
    postfitRMSVals = NaN(2,1);
    stateErrors = NaN(L,N);
    
    %% UKF Initialization Stuff
    
    % Params for Weights
    kappa = 3-L;
    %kappa = 1;
    lambda = alpha^2*(L+kappa)-L;
    gamma = sqrt(L+lambda);
    
    % Weights
    W0m = lambda/(L+lambda);
    W0c = W0m + (1-alpha^2+Beta);
    Wi = 1/(2*(L+lambda));
    Wim = [W0m, repmat(Wi,1,2*L)]; % Create vectors for easy vectorization later
    Wic = [W0c, repmat(Wi,1,2*L)];
    
    %WicMat = diag(Wic);
    
    %% Start loop

    Xhat = X0Ref;
    P = P0;
    
    i = 1; % loop counter
    tprev = tVec(1);
    
    % Progress bar setup
    f = waitbar(0,'1','Name','Running UKF with DMC...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    setappdata(f,'canceling',0);

    while i<=N
        
        if getappdata(f,'canceling')
            break
        end
        
        ti = tVec(i);
        y = measurements(i,3:4)'; % current measurement
        
        %% DMC Noise Matrix
        
        DMCt0 = tprev;
        
        Q = zeros(L,L);
        
        if ~((ti-tprev)>3600)
            Qrr = P0(7:9,7:9)*(1/(3*beta^2)*(ti-DMCt0)^3 - 1/(beta^3)*(ti-DMCt0)^2 + ...
                1/(beta^4)*(ti-DMCt0)*(1-2*exp(-beta*(ti-DMCt0))) + 1/(2*beta^5)*(1-exp(-2*beta*(ti-DMCt0))));
            Qrv = P0(7:9,7:9)*(1/(2*beta^2)*(ti-DMCt0)^2 - 1/(beta^3)*(ti-DMCt0) + ...
                1/(beta^3)*exp(-beta*(ti-DMCt0))*(ti-DMCt0) + 1/(beta^4)*(1-exp(-beta*(ti-DMCt0))) - ...
                1/(2*beta^4)*(1-exp(-2*beta*(ti-DMCt0))));
            Qrw = P0(7:9,7:9)*(1/(2*beta^3)*(1-exp(-2*beta*(ti-DMCt0))) - 1/(beta^2)*exp(-beta*(ti-DMCt0))*(ti-DMCt0));
            Qvv = P0(7:9,7:9)*(1/(beta^2)*(ti-DMCt0) - 2/(beta^3)*(1-exp(-beta*(ti-DMCt0))) + ...
                1/(2*beta^3)*(1-exp(-2*beta*(ti-DMCt0))));
            Qvw = P0(7:9,7:9)*(1/(2*beta^2)*(1+exp(-2*beta*(ti-DMCt0))) - 1/(beta^2)*exp(-beta*(ti-DMCt0)));
            Qww = P0(7:9,7:9)*(1/(2*beta^2)*(1-exp(-2*beta*(ti-DMCt0))));

            Q9 = [Qrr Qrv Qrw;
                Qrv Qvv Qvw;
                Qrw Qvw Qww];

            Q(1:9,1:9) = Q9;
        end
        %% Time Update
        
        % Compute sigma points
        chi = real([Xhat, repmat(Xhat,1,L)+gamma*sqrtm(P), repmat(Xhat,1,L)-gamma*sqrtm(P)]);
        
        if i ~= 1 % If not first step, integrate through, otherwise Xbar doesnt change
            
            stateVec = reshape(chi,L*(2*L+1),[]); % Set up ICs for integration
            
            [~, chi] = ode45(@augmentedDynamics, [tprev ti], stateVec, odeOpts); % Integrate 
            chi = reshape(chi(end,:), L, (2*L+1));
            
            % Compute new mean and covar
            Xbar = sum(Wim.*chi, 2);
            %Xbar = chi*Wim';
        else
            Xbar = Xhat;
        end
        
        % Update 
        chiMinusMean = chi - repmat(Xbar,1,(2*L+1));
        Pbar = (Q + chiMinusMean*diag(Wic)*chiMinusMean');

        % Recompute sigma points using new mean and covar
        chi = real([Xbar, repmat(Xbar,1,L)+gamma*sqrtm(Pbar), repmat(Xbar,1,L)-gamma*sqrtm(Pbar)]);
        
        % Predicted measurements
        stationIndex = measurements(i,2);

        YMat = NaN(yLen, 2*L+1);
        for j = 1:(2*L+1)
            YMat(:,j) = (measurementFunction(chi(:,j), stationIndex));
        end

        ybar = sum(Wim.*YMat, 2); % Weighted avg of predicted measurements w each sigma point
        %ybar = YMat*Wim';

        % Compute new covariances
        YMinusMean = YMat - repmat(ybar,1,(2*L+1));

        Pyy = R + YMinusMean*diag(Wic)*YMinusMean'; % innovation
        Pxy = chiMinusMean*diag(Wic)*YMinusMean'; % cross correlation

        %% Measurement Update

        K = (Pxy/(Pyy)); % Kalman gain

        Xhat = real(Xbar + K*(y - ybar)); % Update state

        P = real(Pbar - K*Pyy*K'); % Update covar

        %% Save histories and residuals and stuff
        
        % Histories
        XHist(:,i) = Xhat;
        PHist(:,:,i) = P;
        sigmaHist(:,i) = sqrt(diag(P));

        % State errors (if given true state vals)
        if truthFlag
            stateErrors(:,i) = XTruth(ti/10+1,:)' - Xhat;
        end

        % Residuals
        prefitResiduals(:,i,stationIndex) = y-ybar;
        postfitResiduals(:,i,stationIndex) = y-measurementFunction(Xhat, stationIndex);
        
        %% Set stuff for next loop iteration and update progress bar
        if mod(i,10)==0
            waitbar(i/N,f,sprintf('Step %i of %i',i,N))
        end
        
        i = i+1;
        tprev = ti;
    
    end
    
    delete(f)
    %% Compute RMS 
    
    %postfitRMSVals = (sum(postfitResiduals.^2,2))./N.^(1/2);
    
    postfitRMSValsAccum = zeros(2,1);
    for i = 1:N
        postfitRMSValsAccum = postfitRMSValsAccum + postfitResiduals(:,i).^2;
    end

    postfitRMSVals = (postfitRMSValsAccum.*(1/N)).^(1/2);
   
    if truthFlag
        RMS3DPosAccum = 0;
        RMS3DVelAccum = 0;
        for i = 1:N
            RMS3DPosAccum = RMS3DPosAccum + stateErrors(1:3,i)'*stateErrors(1:3,i);
            RMS3DVelAccum = RMS3DVelAccum + stateErrors(4:6,i)'*stateErrors(4:6,i);
        end

        RMS3DPos = sqrt(RMS3DPosAccum/N);
        RMS3DVel = sqrt(RMS3DVelAccum/N);

        RMS3D = [RMS3DPos; RMS3DVel];
    else
        RMS3D = [];
    end
    %% Outputs
    
    UKFStruct = struct('stateEstimates', XHist, 'filterCovars', PHist, 'filterSigmas', sigmaHist, 'prefits', ...
        prefitResiduals, 'postfits', postfitResiduals, 'postfitRMS', postfitRMSVals, 'RMS3D', RMS3D, 'stateErrors', stateErrors);

end