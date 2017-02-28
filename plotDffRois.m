function plotDffRois(ext)
% This function will load NS Rel map and allow for ROI selection. It will
% then extract the dFF for ROI accross the 6 trials for both NS and
% gratings and plot the traces for each.

if nargin==0
    fprintf('Select data.mat file...\n\n')
    [ext, pathname] = uigetfile('*data.mat','Load data.mat file containing data structure...');
    cd(pathname);
    fprintf('Loading %s...\n\n', ext);
    load (ext, 'data')
    load([pathname ext(1:end-8), 'RelNat']);
    load([pathname ext(1:end-8), 'RelGrat']);
    load([pathname ext(1:end-8), 'NatImages']);
    load([pathname ext(1:end-8), 'Gratings']);
elseif nargin == 1
    fprintf('Loading files...\n\n');
    load([ext(1:end-8), 'RelNat']);
    load([ext(1:end-8), 'RelGrat']);
    load([ext(1:end-8), 'NatImages']);
    load([ext(1:end-8), 'Gratings']);
elseif nargin==3
    load (ext, 'data')
    fprintf('Loading data to plot...\n');
end

%% Input params
saveFlag        =   1;
cropFlags       = [96, 128, 256];   % Sizes indicating files have not been cropped to remove edge effects/
cropSz          =   4;               % Amount to crop from each edge (eg 4 + 1)
cropSt          =   cropSz+1;
xvals = [13,18,30,35,47,52,64]';

%% Plot Rel Maps

% Crop file if necessary
if ismember(size(RelNat,1), cropFlags)
    NatImages = cropMat(NatImages,5);
    Gratings = cropMat(Gratings,5);
    RelNat = cropMat(RelNat,5);
    RelGrat = cropMat(RelGrat,5);
    dff = cropMat(data.dFF,5);
% else
%     dff = data.dFF;
end
%% Normalize maps
medRelNat=median(median(RelNat));
RelNat(RelNat<medRelNat) = NaN;
minvNat = min(min(RelNat));
maxvNat = max(max(RelNat));
RelNat(isnan(RelNat)) = minvNat-((maxvNat-minvNat)/5);
%RelGrat
meanRelGrat=mean(mean(RelGrat));
RelGrat(RelGrat<meanRelGrat) = NaN;
minvGrat = min(min(RelGrat));
maxvGrat = max(max(RelGrat));
RelGrat(isnan(RelGrat)) = minvGrat-((maxvGrat-minvGrat)/5);
RelNatZ = (RelNat - mean(RelNat(:))) ./ std(RelNat(:));
RelGratZ = (RelGrat - mean(RelNat(:))) ./ std(RelNat(:));

%% Select roi

% Plot NS Rel Map and Max Dff to select roi
ddd=[0 0 0;jet(10)];

f1 = figure; hold on;
imagesc( RelNat ); colormap(ddd); axis square; axis off; title('NS Rel');


% Start roi selection and recover positions
h = impoly;
positions = getPosition(h);
[m,n,z] = size(dff); 
mask = poly2mask(positions(:,1),positions(:,2),m,n);
close (f1)
% Extract dff values for NS 
natblock = [];
for i = 1:6
    curNat = NatImages{i};
    numFrames = size(curNat,3);
    for j = 1:numFrames
        curr_frame = curNat(:,:,j);
        dffROI_NS(i,j) = sum(sum(curr_frame(mask(:,:))))/sum(sum(mask(:,:)));
    end
    natblock  = [natblock dffROI_NS(i,:)];
end


% Extract dff values for Gratings 
gratblock = [];
for i = 1:6
    curGrat = Gratings{i};
    numFrames = size(curGrat,3);
    for j = 1:numFrames
        curr_frame = curGrat(:,:,j);
        dffROI_Grat(i,j) = sum(sum(curr_frame(mask(:,:))))/sum(sum(mask(:,:)));
    end
    gratblock = [gratblock dffROI_Grat(i,:)];
end

xs = repmat(xvals,1,2);
ys = repmat([0 7],length(xvals),1);
f2 = figure; hold on
subplot(2,2,1)
imagesc( RelNat ); colormap(ddd); axis square; axis off; title('NS Rel Map');
patch(positions(:,1),positions(:,2),'red','FaceAlpha',0.75);
subplot(2,2,2)
plotRaster(dffROI_NS ,data.newFrameRate);
ylabel('Trials');
xlabel('Time (s)');
for l = 1:length(xs)
    line(xs(l,:),ys(l,:))
end
subplot(2,2,3)
imagesc( RelGrat ); colormap(ddd); axis square; axis off; title('Gratings Rel Map');
patch(positions(:,1),positions(:,2),'red','FaceAlpha',0.75);
subplot(2,2,4)
plotRaster(dffROI_Grat,data.newFrameRate);
ylabel('Trials');
xlabel('Time (s)');
set(f2,...
        'color'         ,'w'                            ,...
        'Name'          ,(data.filename(1:end-9))       ,...
        'NumberTitle'   ,'off'                          );
    annotation('textbox', [0 0.9 1 0.1]                 ,...
        'String', [data.filename(1:end-9)]              ,...
        'EdgeColor', 'none'                             ,...
        'HorizontalAlignment', 'center'                 );
% clearvars -except dffROI_Grat dffROI_NS positions mask
clearvars xs ys
%% Save masks
rois.masks = mask;
rois.dffNS = dffROI_NS;
rois.dffGrat = dffROI_Grat;
rois.positions = positions;
rois.NatAvgdFF = mean(avgNatResp);
rois.GratAvgdFF = mean(avgGratTrace);
rois.Natblock = reshape(natblock,119,24);
rois.avgNatResp = mean(mean(rois.Natblock(1:91,:),2));
rois.GratBlock = reshape(gratblock,56,72);
rois.avgGratResp = mean(mean(rois.GratBlock(1:14,:),2));

if saveFlag
    ext = input('Enter name suffix: ','s');
    save(['roiMasks_stim' ext '.mat'], 'rois');
    saveas(f2,[f2.Name 'stim_roi_' ext '.jpg'],'jpeg');
    saveas(f2,[f2.Name 'stim_roi_' ext '.eps'],'eps2c');
end
clear rois
