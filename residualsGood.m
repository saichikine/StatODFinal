function postFits = residualsGood(states, measurementFunction, measurements, envParams)

    N = length(measurements);
    
    postFits = NaN(2,N,3);

    % Progress bar setup
    f = waitbar(0,'Initializing...','Name','Computing residuals...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    setappdata(f,'canceling',0);
    
    i = 1;
    
    while i <= N
        
        if mod(i,1)==0
%             waitbar(i/N,f,sprintf('%.2f%%',i/N))
            waitbar(i/N,f,sprintf('Working on: Step %i of %i',i,N))
        end
        
        if getappdata(f,'canceling')
            break
        end

        y = measurements(i,3:4)';
        stationIndex = measurements(i,2);

        yPredicted = measurementFunction(states(i,:)', stationIndex);

        postFits(:,i,stationIndex) = y - yPredicted;
        
        i = i+1;

    end
    
    delete(f)
end