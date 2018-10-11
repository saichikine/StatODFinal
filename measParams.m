function YPredicted = measParams(X, stationIndex)
    
    % Unpack stuff from state
    R = X(1:3); % s/c position
    Rdot = X(4:6); % s/c velocity
    Rs = X(7+(stationIndex-1)*3:7+(stationIndex-1)*3+2); % station position
    if stationIndex==1
        br = X(16); % station range bias
        brdot = X(17); % station range-rate bias
    elseif stationIndex==2
        br = X(18);
        brdot = X(19);
    elseif stationIndex==3
        br = X(20);
        brdot = X(21);
    end
    
    omegaVec = [0; 0; X(24)]; % earth spin rate

    YPredicted = zeros(2,1); 
    
    Rsdot = cross3d(omegaVec, Rs);
    
    rho = sqrt((R - Rs)'*(R - Rs)) + br;
    rhodot = (R - Rs)'*(Rdot - Rsdot)/rho + brdot;
    
    YPredicted(1) = rho;
    YPredicted(2) = rhodot;
end