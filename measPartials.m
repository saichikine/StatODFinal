function [dGdX] = measPartials(X, stationState)

    dGdX = zeros(2,7);

    R = X(1:3); %column vector
    Rs = stationState(1:3); %column vector
    Rdot = X(4:6); %column vector
    Rsdot = stationState(4:6); %column vector
    
    rho = sqrt((R - Rs)'*(R - Rs));
    rhodot = (R - Rs)'*(Rdot - Rsdot)/rho;
    
    dGdX(1,1:3) = (R - Rs)'/rho; %d
    
    dGdX(2,1:3) = (Rdot - Rsdot)'/rho - rhodot/rho*dGdX(1,1:3);
    dGdX(2,4:6) = dGdX(1,1:3);
    
    
end