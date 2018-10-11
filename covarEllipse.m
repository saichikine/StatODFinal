function points = covarEllipse(C)

    n=100; % Number of points around ellipse
    p=0:pi/n:2*pi; % angles around a circle

    [eigvec,eigval] = eig(C); % Compute eigen-stuff
    xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
    x = xy(:,1);
    y = xy(:,2);
    
    points = [x y];

end