% Sai Chikine
% 6080 Project 1
% Classical Kalman Filter - Potter Algorithm

% TO FIX: get rid of 'stationStates' in args - should be more robust than
% this

% REQUIRED INPUTS: 
% 'X0_ref': initial reference trajectory
% 'P0_apriori': a-priori covariance
% 'deltaX0_apriori': a-priori deviation 
% 'R': measurement noise matrix R
% 'measurements': measurement data 
% 'stationStates': station state data
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

function CKFStruct = CKFPotterIterated(X0_ref, P0_apriori, deltaX0_apriori, R, measurements, stationStates, dynamicsFunction, measurementPartialsFunction, measurementFunction, params, varargin)
    
    % Handle optional arguments
    numVarArgs = length(varargin);
    if numVarArgs > 3
        error('CKF:TooManyInputs', 'requires at most 3 optional inputs');
    end
    
    optArgs = {odeset('RelTol', 1e-13, 'AbsTol', 1e-20), 1, 0};
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
    
    % 
    measNum = 2; % assume range and range rate for now

    % Start filter
    iterationCounter = 0;
    
    while iterationCounter < iterationNum
        
        %initialization stuff (Step i=1)
        
        XStar_i = X0Ref; % set X*(ti-1) to X0Ref before inner loop starts
        
        tiprev = 0;
        XStar_iPrev = X0Ref;
        x_ihatPrev = initialGuess;
        PPrev = P0_apriori;
        STM_i = eye(XDim);
        
        % Potter stuff
        WPrev = chol(PPrev);
        
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
            Wibar = STM_i*WPrev;
            
            % try to account for multiple station sightings
            
            numMeasurements = 1; % default setting for number of measurements
            
            % count number of measurements with this time stamp (should be
            % in order)
            if i ~= length(y)
                while tVec(i+1) == tVec(i)
                    numMeasurements = numMeasurements + 1;
                end
            end
            
            for j = 1:numMeasurements
                
                stationIndex = y(i+j-1,2);
                
                for k = 1:measNum % only dealing with one single measurement value
                    
                    measPartialResults = measurementPartialsFunction(XStar_i, stationStates(:, stationIndex, i+j-1));
                    yResults = y(i+numMeasurements-1,3:4)';
                    ymResults = measurementFunction(XStar_i, stationStates(:, stationIndex, i+j-1));
                    
                    Yi = yResults(k);
                
                    HTilde_i = [measPartialResults(k,:)];
                    Yi_m = ymResults(k);
                    Rk = R(k,k);
                    yi = Yi - Yi_m;
                    
                    % Start actual Potter stuff
                    FTilde_i = Wibar'*HTilde_i';
                    ai = inv(FTilde_i'*FTilde_i + Rk);
                    Ki = ai*Wibar*FTilde_i;
                    
                    x_ihat = x_ibar + Ki*(yi - HTilde_i*x_ibar);
                    gamma_i = 1/(1+sqrt(ai*Rk^2));
                    W_i = Wibar - gamma_i*Ki*FTilde_i';
                    P_i = W_i*W_i';
                    
                    x_ibar = x_ihat;
                    Wibar = W_i;
                end
                
            end
            
            % Prep for next iteration
            tiprev = ti;
            XStar_iPrev = XStar_i;
            x_ihatPrev = x_ihat;
            WPrev = W_i;
            
            % Save histories (except for XHist which was taken care of
            % already)
            xHist(:,i) = x_ihat;
            PHist(:,:,i) = P_i;
            sigmaHist(:,i) = sqrt(diag(P_i));
            STMHist(:,:,i) = STM_i;
            
            % Save residuals
            prefitResiduals(:,i) = yi;
            postfitResiduals(:,i) = yi - HTilde_i*x_ihat;
            
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
            postfitRMSValsTotal(:,:,iterationCounter+1) = postfitRMSVals;
        end
        
        fprintf('Just finished iteration %i\n', iterationCounter);
        iterationCounter = iterationCounter + 1;
        
    end
    
    CKFStruct = struct('X_Histories', XHistTotal, 'x_histories', xHistTotal, 'Covar_Histories', PHistTotal, 'Stdevs_Histories', sigmaHistTotal, ...
        'STM_Histories', STMHistTotal, 'Prefit_Residuals', prefitResidualsTotal, 'Postfit_Residuals', postfitResidualsTotal, 'Postfit_RMS_Vals', ...
        postfitRMSValsTotal, 'RMS3D', RMS3D, 'State_Errors', stateErrors, 'Iteration_Counter', iterationCounter);
end

