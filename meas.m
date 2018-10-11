function YPredicted = meas(X, stationState)

    YPredicted = zeros(2,1); 
    
    R = X(1:3);
    Rs = stationState(1:3);
    Rdot = X(4:6);
    Rsdot = stationState(4:6);
    
    rho = sqrt((R - Rs)'*(R - Rs));
    rhodot = (R - Rs)'*(Rdot - Rsdot)/rho;
    
    YPredicted(1) = rho;
    YPredicted(2) = rhodot;
end