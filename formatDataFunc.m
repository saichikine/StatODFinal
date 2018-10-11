function [] = formatDataFunc(inputDataName, outputDataString)
    
    S = load(inputDataName);
    fieldNames = fieldnames(S);
    fieldName = fieldNames{1};
    Obs = S.(fieldName);

    station1Vec = [Obs(:,2), Obs(:,5)];
    station2Vec = [Obs(:,3), Obs(:,6)];
    station3Vec = [Obs(:,4), Obs(:,7)];

    observations = [];

    for i = 1:length(station1Vec)
        if ~isnan(station1Vec(i,1))
            observations = [observations; Obs(i,1), 1, Obs(i,2), Obs(i,5)];
        elseif ~isnan(station2Vec(i,1))
            observations = [observations; Obs(i,1), 2, Obs(i,3), Obs(i,6)];
        elseif ~isnan(station3Vec(i,1))
            observations = [observations; Obs(i,1), 3, Obs(i,4), Obs(i,7)];
        end
    end

    measurements = sortrows(observations,1);

    save(outputDataString, 'measurements');
    
end