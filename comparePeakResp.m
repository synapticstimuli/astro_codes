function comparePeakResp(expNum, rois, stimType)

% Input argments: 
%       expNum is a string identifying the two files to be analyzed together
%       rois indicates the specific rois to analyze, 
%       stim type: 'nm' natural movie, 'osi' gratings response. 


switch nargin
    case 0
        currExp = '005';
        sigROIs = [];
        stimType = 'osi';
    case 1
        currExp = expNum;
        sigROIs = [];
        stimType = 'osi';
    case 2
        currExp = expNum;
        sigROIs = rois;
        stimType = 'osi';
    case 3
        currExp = expNum;
        sigROIs = rois;
        stimType = stimType;
end

%% Load data..............

files = dir(['*' currExp '*rois.mat']);

if numel(files)> 2
    fprintf('Too many files found.\n');
    files = uigetfile('*rois.mat','Select rois.mat files for each channel...',...
        'MultiSelect','on');
    prefname = files{1};
    postfname = files{2};
else
    prefname = files(1).name;
    postfname = files(2).name;
end

preFile = load(prefname,'rois');
postFile = load(postfname,'rois');
load(prefname,'Params');

%% Partition and process data...............
if isempty(sigROIs)
    numROIs = size(preFile.rois.dFFNew,1);
    roiID = [1:numROIs]';
else
    numROIs = numel(sigROIs);
    % Create table row headers
    roiID = sigROIs';
end;
 
% Slect stimuls dependet data...
switch stimType
    case 'osi'
    preONResp = preFile.rois.ONResponse;
    preOFFResp = preFile.rois.OFFResponse;
    postONResp = postFile.rois.ONResponse;
    postOFFResp = postFile.rois.OFFResponse;
    
    case 'nm'
        durSec_NM       = Params.durStim_ns + Params.durBlank_nm;                                   % Duration (s) of NS (ON + OFF)
        onLength        = Params.durStim_ns * preFile.rois.newFrameRate;
        offLength       = Params.durBlank_nm * preFile.rois.newFrameRate;
        durFrames_NM    = durSec_NM * preFile.rois.newFrameRate;
        
        for i = 1:numROIs
            precur = preFile.rois.NatMovies(i,:);
            postcur = postFile.rois.NatMovies(i,:);
            
            prefoo = reshape(precur,durFrames_NM,[]);
            postfoo = reshape(postcur,durFrames_NM,[]);
            
            preON = prefoo(1:offLength,:);
            preOFF = prefoo(offLength+1:end,:);
            postON = postfoo(1:offLength,:);
            postOFF = postfoo(offLength+1:end,:);

            preONResp{i} = preON';
            preOFFResp{i} = preOFF';
            postONResp{i} = postON';
            postOFFResp{i} = postOFF';
        end        
end

% Build table variables...
for i = 1:numROIs
  % Process pre values  
  maxPre(i) = max((mean(preONResp{i},2)));                                  %preFile.rois.maxOriResp(i);
  meanOnPre(i) = mean((mean(preONResp{i},2)));                              %preFile.rois.meanOriResp(i);
  meanOffPre(i) = mean((mean(preOFFResp{i},2)));
  normPre{i} = (mean(preONResp{i},2)-(mean(preOFFResp{i},2)))./(mean(preOFFResp{i},2));
  normPreMean(i) = mean(normPre{i});
  normPreMax(i) = max(normPre{i});
  
  
  % Process post values
  maxPost(i) = max((mean(postONResp{i},2)));                                  %preFile.rois.maxOriResp(i);
  meanOnPost(i) = mean((mean(postONResp{i},2)));                              %preFile.rois.meanOriResp(i);
  meanOffPost(i) = mean((mean(postOFFResp{i},2)));
  normPost{i} = (mean(postONResp{i},2)-(mean(postOFFResp{i},2)))./(mean(postOFFResp{i},2));
  normPostMean(i) = mean(normPost{i});
  normPostMax(i) = max(normPost{i});

end

%% Build table with data......................
MaxResponse = [maxPre' maxPost'];
MeanONResponse = [meanOnPre' meanOnPost']; 
MeanOFFResponse = [meanOffPre' meanOffPost'];
NormMeanResponse = [normPreMean' normPostMean']; 
NormMaxResponse = [normPreMax' normPostMax'];
  
T = [];
T = table;
T.Rois = [num2str(roiID); 'µ:'];
T.MaxResponse = [MaxResponse; mean(MaxResponse)];
T.MeanONResponse = [MeanONResponse; mean(MeanONResponse)];
T.MeanOFFResponse = [MeanOFFResponse; mean(MeanOFFResponse)];
T.NormMeanResponse = [NormMeanResponse; mean(NormMeanResponse)];
T.NormMaxResponse = [NormMaxResponse; mean(NormMaxResponse)];
T.Properties.UserData = {'Experiment files:'; ['preFile: ' prefname]; ['postFile: ' postfname]};
clc;
disp(T.Properties.UserData);
disp(T);
tableName = [prefname(18:end-4), '-table-', stimType, '.txt'];
writetable(T, tableName);

% zscore table
Tz = [];
Tz = table;
Tz.Rois = num2str(roiID);
Tz.MaxResponse = zscore(MaxResponse);
Tz.MeanONResponse = zscore(MeanONResponse);
Tz.MeanOFFResponse =  zscore(MeanOFFResponse);
Tz.NormMeanResponse = zscore(NormMeanResponse);
Tz.NormMaxResponse = zscore(NormMaxResponse);
Tz.Properties.UserData = {'Experiment files:'; ['preFile: ' prefname]; ['postFile: ' postfname]};
Tz.Properties.Description = ['Zscore values for ',tableName];

disp(Tz.Properties.UserData);
disp(Tz);
writetable(Tz,[prefname(18:end-4), '-Ztable-', stimType, '.txt']);

%% Plot data (before-after plots)..................

f = figure('color','w'              ,...
    'position', [200 10 750 700]    ,...
    'defaultaxesfontsize', 14       );
for i = 2:6
    subplot_tight(3,2,i-1,[0.07 0.05]);
    preData = eval(['T.',T.Properties.VariableNames{i},'(:,1)']);
    postData = eval(['T.',T.Properties.VariableNames{i},'(:,2)']);
    h = plotBeforeAfter(preData,postData,{'Pre', '2 wk'});
    title([T.Properties.VariableNames{i}]);
    
end
suptitle(tableName(1:end-4));
saveName = [prefname(18:end-4), '-pairedPlot-', stimType];
saveas(f, saveName, 'eps2c');

%% Plot histograms.....................
numPlots = length(T.Properties.VariableNames);
figure('color','w','Position',[100 100 1000 720]); hold on;
ct = 0;
for i = 2:numPlots
    ct = ct +1;
    subplot_tight(2,numPlots-1,ct,[0.1,0.02]);
    hold on;
    h(1) = histogram(eval(['T.',T.Properties.VariableNames{i},'(:,1)']));
    h(2) = histogram(eval(['T.',T.Properties.VariableNames{i},'(:,2)']));
    title([T.Properties.VariableNames{i}]);
    subplot_tight(2,numPlots-1,ct+(numPlots-1),[0.1,0.02]);
    hold on;
    h(3) = histogram(eval(['Tz.',Tz.Properties.VariableNames{i},'(:,1)']));
    h(4) = histogram(eval(['Tz.',Tz.Properties.VariableNames{i},'(:,2)']));
    title([Tz.Properties.VariableNames{i},' (zscore)']);
    for j = 1:4
        h(j).NumBins = 6;
%         h(j).BinWidth = ;
%         h(j).Normalization = 'probability';
    end
end
suptitle(tableName(1:end-4));

%% For later

% normpostMaxPO(i) = normPre(prePOIdx(i));
  
  % oris = 0:30:180;
  % for j = 1:numel(oris)-1
  %     curPO = preFile.rois.PO(sigROIs(i));
  %     low = oris(j);
  %     upper = oris(j+1);
  %     curPOIdx(j) = curPO>low & curPO<upper;
  % end
  % prePOIdx(i) = oris(curPOIdx);
  
