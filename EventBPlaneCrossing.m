function [position, isterminal, direction] = EventBPlaneCrossing(t, X, SHat)
    position = dot(X(1:3), SHat);
    isterminal = 1;
    direction = 0;
end