function batchStruct = batchFilter(initStuff, measurements, dynamicsFunction, measurementFunction, ...
    measurementPartialsFunction, batchParams, envParams)


    %% Input Conditioning
    x0bar = initStuff{1};
    X0Ref = initStuff{2};
    P0 = initStuff{3};
    R = initStuff{4};
    
    iterNum = batchParams{1};

    L = length(X0Ref); % length of state vector
    tVec = measurements(:,1);
    
    %% Create storage stuff
    
    preFits = NaN(2,length(measurements),3,iterNum);
    postFits = NaN(2,length(measurements),3,iterNum);
    X0Hist = NaN(L,iterNum);
    STMHist = NaN(L,L,length(measurements),iterNum);
    PHist = NaN(L,L,iterNum);

    %% Start iteration loop
    
    iterationCounter = 1;
    
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-18);
    
    while iterationCounter<iterNum

        %initialization stuff (STEP i=1)

        STMiprev = eye(L); %reset STM to identity before inner loop starts
        STMiprev = STMiprev(:);
        
        fprintf('Beginning batch iteration %i... integrating reference trajectory...',iterationCounter)
        % Integrate entire reference trajectory
        ICs = [X0Ref; STMiprev];
        [~,XRef] = ode45(@(t,X) dynamicsFunction(t,X,envParams), tVec, ICs, opts);
        
        fprintf('...done. Beginning processing measurements.\n')
        Lambda = inv(P0); %reset to a priori covar
        if iterationCounter == 1
            N = zeros(25,1);
        else
            N = (P0)\x0bar;
        end
        
        %% Start main loop and create progress bar

        f = waitbar(0,'Starting batch iteration...','Name','Running batch...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

        setappdata(f,'canceling',0);
        
        i = 1;
        while i <= length(measurements)
            
            if mod(i,1)==0
%             waitbar(i/N,f,sprintf('%.2f%%',i/N))
                waitbar(i/length(measurements),f,sprintf('Working on: step %i of %i',i,length(measurements)))
            end
        
            if getappdata(f,'canceling')
                break
            end

            Yi = measurements(i,3:4)';
            XRefState = XRef(i,:)';
            STMi = reshape(XRefState(L+1:end),L,[]);
            
            % Save STM
            STMHist(:,:,i,iterationCounter) = STMi;
            
            % Do measurement stuff
            stationIndex = measurements(i,2);

            Htildei = measurementPartialsFunction(XRefState, stationIndex);
            Yi_m = measurementFunction(XRefState, stationIndex);

            %4th box
            yi = Yi - Yi_m;
            preFits(:,i,stationIndex,iterationCounter) = yi;
            Hi = Htildei*STMi; 
            Lambda = Lambda + Hi'/(R)*Hi;
            N = N + Hi'/(R)*yi;
            
            i = i+1;
            
        end

        x0hat = Lambda\N; %solve normal equations

        X0Ref = X0Ref + x0hat; %Update new reference trajectory
        x0bar = x0bar - x0hat; %Shift new x0bar to be w.r.t. new ref trajectory
        
        % Save
        X0Hist(:,iterationCounter) = X0Ref;
        PHist(:,:,iterationCounter) = inv(Lambda);

        iterationCounter = iterationCounter+1;

    end
    
    batchStruct = struct('stateEstimates', X0Hist, 'filterCovars', PHist, 'filterSigmas', sigmaHist, 'prefits', ...
        prefitResiduals, 'postfits', postfitResiduals, 'postfitRMS', postfitRMSVals, 'RMS3D', RMS3D, 'stateErrors', stateErrors);
end