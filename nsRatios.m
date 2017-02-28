% function nsRatios

fprintf('Select data for first condition...\n')
preDir = uigetdir('Select data for first condition...\n');
fprintf('Select data for second condition...\n')
postDir = uigetdir('Select data for second condition...\n');

% Load first condition data
load([preDir filesep 'data' filesep 'NatRelData']);
load([preDir filesep 'data' filesep 'GratRelData']);
preNatData.slopes = natSlopes;
preGratData.slopes = gratSlopes;
preNatData.aocNat = aocNat;
preGratData.aocGrat = aocGrat;


clearvars -except preDir postDir preNatData preGratData

% Load first condition data
load([postDir filesep 'data' filesep 'NatRelData']);
load([postDir filesep 'data' filesep 'GratRelData']);
postNatData.slopes = natSlopes;
postGratData.slopes = gratSlopes;
postNatData.aocNat = aocNat;
postGratData.aocGrat = aocGrat;

clearvars -except preDir postDir preNatData preGratData postNatData postGratData

dataL = length(preNatData.slopes);

for i = 1:dataL
    preSlopes(i) = preNatData.slopes(i) / preGratData.slopes(i);
    postSlopes(i) = postNatData.slopes(i) / postGratData.slopes(i);
    preAreas(i) = preNatData.aocNat(i) / preGratData.aocGrat(i);
    postAreas(i) = postNatData.aocNat(i) / postGratData.aocGrat(i);
end

%% Plot ratios

f = figure; hold on;
subplot(1,2,1)
boxScatter([abs(preSlopes') abs(postSlopes')], f);
title ('Slope Ratios');
subplot(1,2,2)
boxScatter([abs(preAreas') abs(postAreas')], f);
title('Area Ratios');
