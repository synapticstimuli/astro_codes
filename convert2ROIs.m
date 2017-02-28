function convert2ROIs(dataFile)

%% Load files.........

switch nargin
    case 0
        [dataFile, pathname] = uigetfile('*data.mat', 'Select data.mat file to convert...');
        if dataFile == 0
            fprintf('No file selected. Script aborted...\n');
            return
        end
        load(dataFile);
    case 1
        fprintf('Loading data structure...\n');
        load(dataFile);
        pathname = [pwd filesep];
end


%% Convert data structure to rois.....

% Copy basic info
rois = [];
rois.WindowSize = data.WindowSize;
rois.stimProtocol = data.stimProtocol;
rois.stimDur = data.stimDur;
rois.windowSize = data.WindowSize;
rois.tifName = data.filename;
rois.numFrames = data.numFrames;
rois.xPixels = data.xPixels;
rois.yPixels = data.yPixels;
rois.rawFrameRate = data.rawFrameRate;
rois.newFrameRate = data.newFrameRate;

% Extract grid dff structure data in to roi based structure....
NumBins = rois.xPixels / rois.windowSize;
ct = 0;
tic;
rois.dFFNew = zeros(NumBins*NumBins,length(data.dFF));
for i = 1:NumBins
    for j = 1:NumBins
        ct = ct+1;
        clc
        fprintf('Converting DFF...\n');
        fprintf('Completed: %0.2f%% (%0.2fs)\n ',(ct/(NumBins*NumBins))*100,toc);
        rois.dFFNew(ct,:) = data.dFF(i,j,:);
    end
end

% Convert Reliability data
try
    relNMfile = strrep(dataFile,'_data.mat','_RelNat.mat');
    relGratfile = strrep(dataFile,'_data.mat','_RelGrat.mat');
    fprintf('Loading reliability files...\n');
    load(relNMfile);
    load(relGratfile);
catch
    relNMfile = uigetfile('*RelNat.mat','Select NM rel file...');
    relGratfile = uigetfile('*RelGrat.mat','Select NM rel file...');
    load(relNMfile);
    load(relGratfile);
end
rois.Reliability = reshape(RelNat',1,[]);
rois.ReliabilityGrat = reshape(RelGrat',1,[]);  

% Create roi positon matrix
fprintf('Creating roi position matrix...\n\n');
rois.positions = makeGridMask([rois.xPixels, rois.yPixels], rois.windowSize);

%% Save.................

% Save new file and variables................
fprintf('Saving...\n\n');
saveName = strrep(dataFile,'_data.mat','-rois.mat');
save([pathname saveName], 'rois', '-v7.3');

%% Clear and close all
% ca;