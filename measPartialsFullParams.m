function [dGdX] = measPartialsFullParams(X, stationIndex)

    dGdX = zeros(2,25);

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
    
    omegaVec = [0; 0; X(24)];
    Rsdot = cross3d(omegaVec,Rs);
    
    rho = sqrt((R - Rs)'*(R - Rs)) + br;
    rhodot = (R - Rs)'*(Rdot - Rsdot)/rho + brdot;
    
    dGdX(1,1:3) = (R - Rs)'/rho; %d
    
    dGdX(2,1:3) = (Rdot - Rsdot)'/rho - rhodot/rho*dGdX(1,1:3);
    dGdX(2,4:6) = dGdX(1,1:3);
    
    if stationIndex == 1
        dGdX(:,7:9) = -dGdX(:,1:3);
        dGdX(:,16:17) = ones(2,2);
    elseif stationIndex == 2
        dGdX(:,10:12) = -dGdX(:,1:3);
        dGdX(:,18:19) = ones(2,2);
    elseif stationIndex == 3
        dGdX(:,13:15) = -dGdX(:,1:3);
        dGdX(:,20:21) = ones(2,2);
    end
    
end