function [position, isterminal, direction] = EventRSOI3(t, X, RSOI3)
    position = norm(X(1:3))-RSOI3;
    isterminal = 1;
    direction = -1;
end