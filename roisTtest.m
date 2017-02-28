% Run ttest for on and off periods for NM and Gratings

if ~exist('rois','var')
    error('Nor rois var foud...');
end

NumROIs = numel(rois.Label);

% Stimulus-specific trial structures...
durSec_NM       = Params.durStim_ns + Params.durBlank_nm;                                   % Duration (s) of NS (ON + OFF)
durFrames_NM    = durSec_NM * rois.newFrameRate;                                       % Num frames of NS stim (ON + OFF)
durTrials_NM    = durFrames_NM * Params.repMovies;                                     % Total length of NS (ON + OFF) per trial
durSec_Grats    = Params.durStim_grats + Params.durBlank_grats;                        % Duration (s) of Grats (ON + OFF)
durFrames_Grats = durSec_Grats * rois.newFrameRate;                                 % Num frames of Grats stim (ON + OFF)
durTrials_Grats = durFrames_Grats * Params.numOris;
numTrials_Grats = Params.numOris * Params.repTrials;  
%% Run ttest for NM...
clc;
fprintf('Running ttest for NM...\n');
progressbar;

% Calculate p value for responses to NM
for i = 1:NumROIs
    cur = rois.NatMovies(i,:);
    foo = reshape(cur,durFrames_NM,[]);
    ON = foo(1:20,:);
    OFF = foo(21:end,:);
    [h, p] = ttest2(ON(:), OFF(:));
    rois.sigNMResp(i) = p;
    progressbar(i/NumROIs);
    h_pb = gcf;
end
close(h_pb);

%% Run ttest for gratings...
progressbar;
fprintf('Running ttest for gratings...\n');

for i =  1:NumROIs
    curON = rois.ONResponse{i}(:);
    curOFF = rois.ONResponse{i}(:);
    [h p] = ttest2(curON,curOFF);
    rois.sigOriResp(i) = p;
    progressbar(i/NumROIs);
    h_pb = gcf;
end

close(h_pb);

%% Save
fprintf('Saving...\n');
save(rois.fname,'rois','-append');