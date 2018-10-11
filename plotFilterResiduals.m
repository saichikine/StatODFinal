function [] = plotFilterResiduals(measurements, prefits, postfits, startStopIndices, varargin);

    %% Handle optional arguments
    numVarArgs = length(varargin);
    if numVarArgs > 2
        error('plotResiduals:TooManyInputs', 'requires at most 2 optional arguments');
    end
    
    optArgs = {[], []};
    optArgs(1:numVarArgs) = varargin;
    [ylimRange, ylimRangeRate] = optArgs{:};
    %%
    startIndex = startStopIndices(1);
    stopIndex = startStopIndices(2);
    
    figure
    subplot(211); hold on
    r1 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(1,:,1), '.','DisplayName', 'Range Pre-fits, Station 1');
    r2 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(1,:,2), '.','DisplayName', 'Range Pre-fits, Station 2');
    r3 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(1,:,3), '.','DisplayName', 'Range Pre-fits, Station 3');
    b1 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
    b2 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [-3*5/1000, -3*5/1000], 'k--'); hold off
    if nargin>4
        ylim(ylimRange)
    end
    legend([r1 r2 r3 b1])
    ylabel('[km]')
    title('Pre-fits');
    subplot(212); hold on
    rr1 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(2,:,1), '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
    rr2 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(2,:,2), '.','DisplayName', 'Range-Rate Pre-fits, Station 2');
    rr3 = scatter(measurements(startIndex:stopIndex,1)./86400, prefits(2,:,3), '.','DisplayName', 'Range-Rate Pre-fits, Station 3');
    b1 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
    b2 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
    if nargin>4
        ylim(ylimRangeRate)
    end
    legend([rr1 rr2 rr3 b1])
    ylabel('[km/s]')
    xlabel('Time [days]')

    figure
    subplot(211); hold on
    r1 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(1,:,1), '.','DisplayName', 'Range Pre-fits, Station 1');
    r2 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(1,:,2), '.','DisplayName', 'Range Pre-fits, Station 2');
    r3 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(1,:,3), '.','DisplayName', 'Range Pre-fits, Station 3');
    b1 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
    b2 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [-3*5/1000, -3*5/1000], 'k--'); hold off
    if nargin>4
    	ylim(ylimRange)
    end
    legend([r1 r2 r3 b1])
    ylabel('[km]')
    title('Post-fits');
    subplot(212); hold on
    rr1 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(2,:,1), '.','DisplayName', 'Range-Rate Post-fits, Station 1');
    rr2 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(2,:,2), '.','DisplayName', 'Range-Rate Post-fits, Station 2');
    rr3 = scatter(measurements(startIndex:stopIndex,1)./86400, postfits(2,:,3), '.','DisplayName', 'Range-Rate Post-fits, Station 3');
    b1 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
    b2 = plot([measurements(startIndex,1), measurements(stopIndex,1)]./86400, [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
    if nargin>4
        ylim(ylimRangeRate)
    end
    legend([rr1 rr2 rr3 b1])
    ylabel('[km/s]')
    xlabel('Time [days]')
    
end