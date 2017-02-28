function [ch1Mat, ch2Mat] = roiOverlap(varargin)
filtPscore = 0;
saveFlag = 1;
discAngles = 30;
edges = 0:discAngles:359;
angThr = 30;
alphaVal = .3;

%% Get files
switch numel(varargin)
    case 0
        testVar = input('Select the response feature to base the overlap map on (OSI, Rel, PO):\n','s');
    case 1
        testVar = varargin{1};
    case 2
        testVar = varargin{1};
        curDir = varargin{2};
        cd(curDir);
end

%%
files = dir('*rois.mat');
if numel(files)> 2
    fprintf('Too many reg2.tif files found.\n');
    files = uigetfile('*rois.mat','Select rois.mat files for each channel...',...
        'MultiSelect','on');
    ch1fname = files{1};
    ch2fname = files{2};
else
    ch1fname = files(1).name;
    ch2fname = files(2).name;
end


%% Load roi files...............................
fprintf('Loading files...');
ch1 = load(ch1fname,'rois');
ch2 = load(ch2fname,'rois');
fprintf('done!\n');

%% Load values for each channel............................................
switch testVar
    case 'OSI'
        ch1vals = ch1.rois.OSI;
        ch2vals = ch2.rois.OSI;
        ch1pvals = ch1.rois.sigOriResp;
        ch2pvals = ch2.rois.sigOriResp;
        suptitleTxt = 'OSI';
        testCase = 0;
    case 'Rel'
        ch1vals = ch1.rois.Reliability;
        ch2vals = ch2.rois.Reliability;
        ch1pvals = ch1.rois.sigNMResp;
        ch2pvals = ch2.rois.sigNMResp;
        suptitleTxt = 'NM-Rel';
        testCase = 0;
    case 'PO'
        ch1vals = ch1.rois.PO;
        ch2vals = ch2.rois.PO;
        
%         ch1vals(ch1pvals>0.05) = NaN;
%         ch2vals = ch2.rois.PO(ch2pvals<=0.05);
        suptitleTxt = 'Preferred Orientation';
        testCase = 1;
end
%% Create colocalization matrix
matSize = ch1.rois.xPixels / ch1.rois.windowSize;
expName = ch1fname(1:strfind(ch1fname,'-Ch1')-1);

switch testCase
    case 0      % Process OSI or Reliability values for overlap....
        if ~filtPscore
            maxVal = 2;
            fprintf('Using vals with zscore > %1.3f...\n',maxVal);
            ch1Idx = zscore(ch1vals)>maxVal;
            ch2Idx = zscore(ch2vals)>maxVal;
            testTxt = ['(zscore)'];
        elseif filtPscore
            maxVal = 0.05;
            fprintf('Using pvals < %1.3f...\n', maxVal);
            ch1Idx = ch1pvals<maxVal;
            ch2Idx = ch2pvals<maxVal;
            testTxt = ['(p', num2str(maxVal), ')'];
        end
    case 1      % Process PO values for overlap....
        values = [];
        for i = 2: length(edges)
            values{i-1} = num2str(edges(i));
        end
        ch1Idx = str2num(char(discretize(ch1vals,edges,'categorical',values)'));
        ch2Idx = str2num(char(discretize(ch2vals,edges,'categorical',values)'));
        simIdx = abs(ch1Idx - ch2Idx) == 0;
        ch1Idx = (ch1Idx/180);%.* simIdx;
        ch2Idx = (ch2Idx/180);%.* simIdx;
        ch1Idx(ch1pvals>0.05) = 0;
        ch2Idx(ch2pvals>0.05) = 0;
%         poMap = cmap_angles(180, discAngles);
end


    ch1Mat = reshape(ch1Idx,matSize,[])';
    ch2Mat = reshape(ch2Idx,matSize,[])';
% Create rgb matrix with thresholded channel data
    rgbImage = zeros(matSize, matSize, 3);
    rgbImage(:,:,1) = ch1Mat;
    rgbImage(:,:,2) = ch2Mat;

%% Plot figure...........................................................
f = figure('color','w');
% rgbImage = alphaVal * rgbImage;
imagesc(rgbImage);
axis square; box off;
title(expName);
suptitle([suptitleTxt, '-', testTxt]);
set(gca,                            ...
    'XMinorGrid','on',...
    'YMinorGrid','on',...
    'GridColor',[1 1 1],...
    'GridAlpha',0.5,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'MinorGridColor',[1 1 1],...
    'MinorGridLineStyle','-',...
    'MinorGridAlpha',0.5,...
    'Layer','top'...
    );
% myColorbar(180,6);


%% Save figure

if saveFlag
    saveName = [pwd filesep expName '-roiOverlap-', suptitleTxt, '-', testTxt];
    saveas(gcf, saveName, 'jpeg');
    saveas(gcf, saveName, 'eps2c');
end

%% Nested functions
    function grid = gridlines(matSize)
        grid = zeros(matSize, matSize);
        gLines = 10:10:matSize-1;
        for i = 1:numel(gLines)
            grid(gLines(i),:) = 1;
            for j = 1:numel(gLines)
                grid(:,gLines(i)) = 1;
            end
        end
    end

    function [zVals, zidx] = zscr(inputVals, zThr)
        if nargin == 1
            zThr = 2;
        end
        z = zscore(inputVals);
        zidx = z>zThr;
        zVals = inputVals(zidx);
    end

    function poMat = plotPO(vals)
        maxval = max(vals);
        numVals = numel(vals);
        cc = 255 * colormap('jet');
        roiColor = ceil(length(cc)*(vals ./ maxval));
        desiredColor = zeros(numVals,3);
        poMat = zeros(1,3600,3);
        for i = 1:numVals
            desiredColor(i,:) = cc(roiColor(i),:);
            desiredColor = round(desiredColor);
%             poMat(1,i,:) = desiredcolor(i) ;
        end
        poMat = rshape(poMat,60,60,[]);
        
    end

    function [map] = cmap_angles(maxAngle, interval)
        switch nargin
            case 0
                maxAngle = 180;
                interval = 30;
            case 1
                interval = maxAngle/12;
        end
        map = ([0:interval:maxAngle; 0:interval:maxAngle; 0:interval:maxAngle]')./maxAngle;
    end
end
