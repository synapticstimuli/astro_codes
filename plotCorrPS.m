function plotCorrPS

%% Input params
rVal = 0.5;
ccFolder = 'C:\Users\rodrigog\Desktop\ccData\mecp2-wt';
plotPS = 0;
plotHist = 1;
saveFlag = 1;
%% Load files
% Load vis response data
try
   pathname = [pwd filesep]; 
   fname = dir([pathname filesep '*NatImages.mat']);
   fname =fname.name;
   fprintf('Loading %s:...\n', fname(1:end-23));
   fprintf('NS data...\n');
   load([pathname fname]);
   fprintf('Gratings data...\n');
   load([pathname fname(1:end-13), 'Gratings']);
catch
    fprintf('Select NatImages data.mat file...\n\n')
    [fname, pathname] = uigetfile('*NatImages.mat','Load NatImages.ma file...');
    cd(pathname);
    fprintf('Loading %s:...\n', fname(1:end-23));
    fprintf('NS data...\n');
    load([pathname fname]);
    fprintf('Gratings data...\n');
    load([pathname fname(1:end-13), 'Gratings']);
end
try 
    pathnameSpon = strrep(pathname,'-NS-','-spon-');
    cd(pathnameSpon);
    fnameSpon = dir('*data.mat');
    fnameSpon = fnameSpon.name;
    fprintf('Spontaneous data...\n');
    load (fnameSpon, 'data')
    expName = fnameSpon(1:end-18);
    sponDff = data.dFF;
catch
    % Load spontaneous data
    fprintf('Select spontaneous data.mat file...\n\n')
    [fnameSpon, pathnameSpon] = uigetfile('*data.mat','Load data.mat file containing spontaneous data...');
    cd(pathnameSpon);
    fprintf('Spontaneous data...\n');
    load (fnameSpon, 'data')
    expName = fnameSpon(1:end-18);
    sponDff = data.dFF;
end

    
%% Extract dff data    
% NS data
natdff = [];
for i =1:6
    natdff = cat(3,natdff,NatImages{i});
end
natdff = cropMat(natdff,5);
% natRespDff = repmat(natdff,24,119);
% Grat data
gratdff = [];
for i =1:6
    gratdff = cat(3,gratdff,Gratings{i});
end
gratdff = cropMat(gratdff,5);
% Spontaneous data
sponDff = cropMat(sponDff,5);

%% Run correlation and thresholding

[natCCImg,natThreshImg] = CrossCorrImage(rVal, natdff);
set(gcf, 'Name',[expName '-NatCC']);
[gratCCImg,gratThreshImg] = CrossCorrImage(rVal, gratdff);
set(gcf, 'Name',[expName '-GratCC']);
[sponCCimg, sponTreshImg] = CrossCorrImage(rVal, sponDff);
set(gcf, 'Name',[expName '-SponCC']);

%% Plot PS of cc maps 
if plotPS
    % h1 = figure; hold on
    % plotMultiPS({gratThreshImg,natThreshImg,sponTreshImg});
    dataSets.mats = {gratThreshImg,natThreshImg,sponTreshImg};
    dataSets.titles = {'Grat','NM','Spon'};
    for i = 1:length(dataSets.mats)
        plotMultiPS(dataSets.mats{i})
        title(['PS - ' dataSets.titles{i}]);
        legend('off');
    end
end

%% Histogram of pixels in cc maps
if plotHist
    f2 = figure;
    hst(1) = subplot(1,3,1);
    h2 = histogram(natCCImg(natCCImg>=0.02));
    title('Nat Movies');
    hst(2) = subplot(1,3,2);
    h3 = histogram(gratCCImg(gratCCImg>=0.02));
    title('Gratings');
    hst(3) = subplot(1,3,3);
    h4 = histogram(sponCCimg(sponCCimg>=0.02));
    title('Spontaneous');
    ylims = max(cell2mat(get(hst,'YLim')));
    xlims = max(cell2mat(get(hst,'XLim')));
    for i=1:3;
        set(hst(i),'YLim',ylims,'XLim',xlims);
    end
    set(f2,...
        'color'         ,'w'                ,...
        'Name'          ,[expName '-hist']  ,...
        'Position', [101 544 1600 270]) ;
    
    f3 = figure; %hold on
    hst(1) = subplot(1,3,1);
    h5 = histogram(natCCImg(natCCImg>=rVal));
    title('Nat Movies');
    hst(2) = subplot(1,3,2);
    h6 = histogram(gratCCImg(gratCCImg>=rVal));
    title('Gratings');
    hst(3) = subplot(1,3,3);
    h7 = histogram(sponCCimg(sponCCimg>=rVal));
    title('Spontaneous');
    ylims = max(cell2mat(get(hst,'YLim')));
    xlims = max(cell2mat(get(hst,'XLim')));
    for i=1:3;
        set(hst(i),'YLim',ylims,'XLim',xlims);
    end
    set(f3,...
        'color'         ,'w'                ,...
        'Name'          ,[expName '-histInset']  ,...
        'Position', [101 544 1600 270]) ;
end


%% Save
ccMaps.name = expName;
ccMaps.natCC = natCCImg;
ccMaps.gratCC = gratCCImg;
ccMaps.sponCC = sponCCimg;
ccMaps.natThr = natThreshImg;
ccMaps.gratThr = gratThreshImg;
ccMaps.sponThr = sponTreshImg;
if saveFlag
    cd(ccFolder);
    save([expName '.mat'],'ccMaps');
    fprintf('Saving figures...\n');
    t = get(0,'children');
    for i=1:length(t)
        saveas(t(i),[t(i).Name '.jpg'],'jpeg')
        saveas(t(i),[t(i).Name '.eps'],'eps2c')
    end
end
% pause(3)
% close all
% 
cd(pathname);
cd(['..' filesep])
fprintf('Finished with exp: %s:\n', fname(1:end-23));
end