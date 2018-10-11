% Sai Chikine
% 6080 Project 1
% Classical Kalman Filter with State Noise Compensation

% TO FIX: get rid of 'stationStates' in args - should be more robust than
% this

% REQUIRED INPUTS: 
% 'X0_ref': initial reference trajectory
% 'P0_apriori': a-priori covariance
% 'deltaX0_apriori': a-priori deviation 
% 'R': measurement noise matrix R
% 'measurements': measurement data 
% 'stationStates': station state data
% 'Q': Q matrix for SNC
% 'dt': time step/delta t for SNC
% 'dynamicsFunction': dynamics model (with STM) 
% 'measurePartialsFunction': measurement partials function
% 'measurementFunction': expected measurement function
% 'params': parameters needed
%
% OPTIONAL INPUTS:
% 'odeOpts': odeset options for integrator
% 'iterationNum': number of iterations for batch filter
%
% OUTPUTS:
% 'CKFStruct': struct consisting of:
%   'X_Histories': X history over all time steps and iterations
%   'x_Histories': x (deviation) history over all time steps and iterations
%   'Covar_Histories': P history over all time steps and iterations
%   'Stdevs_Histories': standard deviation (sqrt of diagonal elements of P) history over all time steps and
%   iterations
%   'STM_Histories': STM history over all time steps and iterations
%   'Prefit_Residuals': prefit residuals over all time steps and iterations
%   'Postfit_Residuals': postfit residuals over all time steps and
%   iterations
%   'Postfit_RMS_Vals': componentwise RMS values over all time steps and
%   iterations
%   'Iteration_Counter': number of iterations of CKF
%
% ***Note that this is an iterated batch filter. Defaults to 1 iteration.

function CKFStruct = CKFGeneral(X0_ref, P0_apriori, deltaX0_apriori, R, measurements, ...
    stationStates, Q, dt, dynamicsFunction, measurementPartialsFunction, measurementFunction, params, varargin)
    
    % Handle optional arguments
    numVarArgs = length(varargin);
    if numVarArgs > 3
        error('CKF:TooManyInputs', 'requires at most 3 optional inputs');
    end
    
    optArgs = {odeset('RelTol', 1e-13, 'AbsTol', 1e-13), 1, 0};
    optArgs(1:numVarArgs) = varargin;
    [odeOpts, iterationNum, XTruth] = optArgs{:};
    
    if XTruth == 0
        truthFlag = 0;
    else
        truthFlag = 1;
    end
    
    % Input Conditioning
    XDim = length(X0_ref); %length of state vector
    y = measurements; % measurement data (change name to y for brevity, keep as measurements in args for clarity)
    N = length(y); %total number of measurements
    tVec = y(:,1); % assume first column of meas data is time stamp
    X0Ref = X0_ref; % reference/nominal from input
    initialGuess = deltaX0_apriori; % initial deviation from input

    % Start CKF
    iterationCounter = 0;
    
    while iterationCounter < iterationNum % Number of iterations of batch (default is 1)
        
        %initialization stuff (Step i=1)
        
        XStar_i = X0Ref; % set X*(ti-1) to X0Ref before inner loop starts
        
        tiprev = 0;
        XStar_iPrev = X0Ref;
        x_ihatPrev = initialGuess;
        PPrev = P0_apriori;
        STM_i = eye(XDim);
        
        % Set up stuff to be stored
        XHist = zeros(XDim,N);
        XHist(:,1) = XStar_i;
        xHist = zeros(XDim,N);
        STMHist = zeros(XDim,XDim,N);
        PHist = zeros(XDim,XDim,N);
        sigmaHist = zeros(XDim,N);
        stateErrors = NaN(XDim,N);
        
        % Set up residuals
        prefitResiduals = zeros(2,N);
        postfitResiduals = zeros(2,N);
        postfitRMSVals = zeros(2,N); %componentwise RMS values
        
        % Precompute Gamma*Q*Gamma^T
        
        Gamma = [dt/2*eye(3,3); eye(3,3)];
        %SNCMat = dt^2*Gamma*Q*Gamma';
        SNCMat = zeros(7,7);
        
        i = 1; % loop counter

        while i <= N
            
            ti = tVec(i);
            
            if ti ~= 0 % if at time t0, no need to integrate
                tSpan = [tiprev ti];
                ICs = [XStar_iPrev; reshape(eye(XDim), XDim^2, [])];
                
                [t, XSTM] = ode45(@(t,X)dynamicsFunction(t, X, params), tSpan, ICs, odeOpts);
                STM_i = reshape(XSTM(end,(XDim+1):end), XDim, []); 
                XStar_i = XSTM(end,1:XDim)';
                
                % Save X here (if not, first iteration won't work since
                % first iteration won't have this integration
                XHist(:,i) = XStar_i;
            end
            
            % Time Update
            x_ibar = STM_i*x_ihatPrev;
            
            % Turn off SNC if there's a gap in measurements
            if (ti - tiprev) > 10
                P_ibar = STM_i*PPrev*STM_i';
            else
                P_ibar = STM_i*PPrev*STM_i' + SNCMat; % with SNC (Gamma*Q*Gamma^T)
            end
            
            % try to account for multiple station sightings
            HTilde_i = [];
            Yi = [];
            Yi_m = [];
            R_i = [];
            
            numMeasurements = 1; % default setting for number of measurements
            
            % count number of measurements with this time stamp (should be
            % in order)
            if i ~= length(y)
                while tVec(i+1) == tVec(i)
                    numMeasurements = numMeasurements + 1;
                end
            end
            
            %fprintf('Number of Measurements is %d\n', numMeasurements)
            
            for j = 1:numMeasurements
                if y(i+j-1,2) == 1
                    stationIndex = 1;
                elseif y(i+j-1,2) == 2
                    stationIndex = 2;
                elseif y(i+j-1,2) == 3
                    stationIndex = 3;
                end
                
                HTilde_i = [HTilde_i; measurementPartialsFunction(XStar_i, stationStates(:, stationIndex, (i+j-1)))];
                Yi = [Yi; y(i+numMeasurements-1,3:4)'];
                Yi_m = [Yi_m; measurementFunction(XStar_i, stationStates(:, stationIndex, (i+j-1)))];
                R_i = blkdiag(R_i,R);
            end
            
            yi = Yi - Yi_m;
            Ki = P_ibar*HTilde_i'/(HTilde_i*P_ibar*HTilde_i' + R_i);
            
            % Measurement Update
            x_ihat = x_ibar + Ki*(yi - HTilde_i*x_ibar);
            P_i = (eye(XDim) - Ki*HTilde_i)*P_ibar*(eye(XDim) - Ki*HTilde_i)' + Ki*R_i*Ki'; % Joseph covar update

            % Prep for next iteration
            tiprev = ti;
            XStar_iPrev = XStar_i;
            x_ihatPrev = x_ihat;
            PPrev = P_i;
            
            if ti~=0
                % Save histories (except for XHist which was taken care of
                % already)
                xHist(:,i) = x_ihat;
                PHist(:,:,i) = P_i;
                sigmaHist(:,i) = sqrt(diag(P_i));
                STMHist(:,:,i) = STM_i;

                if truthFlag
                    stateErrors(:,i) = XTruth(ti/10+1,:)' - (XStar_i + x_ihat);
                end

                % Save residuals
                prefitResiduals(:,i) = yi;
                postfitResiduals(:,i) = yi - HTilde_i*x_ihat;
            end
            %fprintf('i is %d\n', i)
            i = i + numMeasurements; %skip to next new time step
        end
        
        % Integrate backwards to t=0 to update new reference trajectory
        CKFFinalEstimate = XHist(1:XDim,end) + xHist(1:XDim,end);
        [~, XOrbitBackwards] = ode45(@(t,X) SRP3BDynamics(t, X, params), [tVec(end), 0], CKFFinalEstimate, odeOpts);
        
        % Update reference trajectory and a priori deviation to be w.r.t.
        % to this new ref 
        initialGuess = initialGuess - (XOrbitBackwards(end,:)' - X0Ref);
        X0Ref = XOrbitBackwards(end,:)';
%         
        % Compute RMS Values
        postfitRMSValsAccum = zeros(2,1);
        for i = 1:N
            postfitRMSValsAccum = postfitRMSValsAccum + postfitResiduals(:,i).^2;
        end
        
        postfitRMSVals = (postfitRMSValsAccum.*(1/N)).^(1/2);
        
        % 3D RMS - Position and Velocity Separately
        
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
        
        % Histories and residuals over iterations
        if iterationCounter == 0
            % If only 1st iteration, don't bother with adding another
            % dimension to the results
            XHistTotal = XHist;
            xHistTotal = xHist;
            PHistTotal = PHist;
            sigmaHistTotal = sigmaHist;
            STMHistTotal = STMHist;
            prefitResidualsTotal = prefitResiduals;
            postfitResidualsTotal = postfitResiduals;
            postfitRMSValsTotal = postfitRMSVals;
            
        elseif iterationCounter ~= 0
            % Grow saved things along another dimension if iterated more
            % than once
            XHistTotal(:,:,iterationCounter+1) = XHist;
            xHistTotal(:,:,iterationCounter+1) = xHist;
            PHistTotal(:,:,:,iterationCounter+1) = PHist;
            sigmaHistTotal(:,:,iterationCounter+1) = sigmaHist;
            STMHistTotal(:,:,:,iterationCounter+1) = STMHist;
            prefitResidualsTotal(:,:,iterationCounter+1) = prefitResiduals;
            postfitResidualsTotal(:,:,iterationCounter+1) = postfitResiduals;
            %postfitRMSValsTotal(:,:,iterationCounter+1) = postfitRMSValsTotal;
        end
        
        iterationCounter = iterationCounter + 1;
        
    end
    
    CKFStruct = struct('X_Histories', XHistTotal, 'x_histories', xHistTotal, 'Covar_Histories', PHistTotal, 'Stdevs_Histories', sigmaHistTotal, ...
        'STM_Histories', STMHistTotal, 'Prefit_Residuals', prefitResidualsTotal, 'Postfit_Residuals', postfitResiduals, 'Postfit_RMS_Vals', ...
        postfitRMSValsTotal, 'RMS3D', RMS3D, 'State_Errors', stateErrors, 'Iteration_Counter', iterationCounter, 'SNCMat', SNCMat);
end

