function plotVizRelMapPS_RG(fname, RelNat, RelGrat)
%% Plots the reliability maps for gratings and ns stims, as well as 
%  the powers spectra analysis. Saves all graphs as eps and jpgs. From RVR.
global dataDir

if nargin==0
    fprintf('Select data.mat file...\n\n')
    [fname, pathname] = uigetfile('*data.mat','Load data.mat file containing data structure...');
    cd(pathname);
    fprintf('Loading %s...\n\n', fname);
    load (fname, 'data')
    load([pathname fname(1:end-8), 'RelNat']);
    load([pathname fname(1:end-8), 'RelGrat']);
    load([pathname fname(1:end-8), 'NatImages']);
    load([pathname fname(1:end-8), 'Gratings']);
elseif nargin==3
    load (fname, 'data')
    fprintf('Loading data to plot...\n');
end

% Variables 
axisMod         = [2 11];            % Range of x axis for plotPowerSpectra, leave empty for range based on data set
saveFlag        =   0;
plotFlag        =   1;
cropFlags       = [96, 128, 256];
cropSz         =   4;
cropSt          =   cropSz+1;          % Amount to crop from each edge
copyFlag        = 0;

%% Check for map sizes, crop if necessary

% If original tiff files were not cropped to remove edge effects, crop the
% map files before plotting and estimating PS values

if ismember(size(RelNat,1), cropFlags)
    fprintf('Cropping movie files...\n');
    % Crop image matrix to remove registration edge effects.
    RelNat = RelNat(cropSt:end-cropSz,cropSt:end-cropSz);
    RelGrat = RelGrat(cropSt:end-cropSz,cropSt:end-cropSz);
    dataDFF = data.dFF(cropSt:end-cropSz,cropSt:end-cropSz,:);
    fprintf('Movie files cropped...\n');
else
    dataDFF = data.dFF;
end

%% Noralize maps

% RelNat
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

%Load activity maps
natActMap=[];
for i=1:length(NatImages)
    natActMap = cat(3,natActMap, NatImages{i});
end
gratNatMap=[];
for i=1:length(Gratings)
    gratNatMap = cat(3,gratNatMap, Gratings{i});
end

normRelNat = RelNat/(vizMaxProjection(natActMap));
normRelGrat = RelGrat/(vizMaxProjection(gratActMap));

%Load avg proj image
try
    projFile = ['AVG_',fname(1:end-9),'.jpg'];
    avgProj = imread(projFile);
catch
    projFile = uigetfile('*.jpg','Select projection jpg file...');
    avgProj = imread(projFile);
end

h1=figure; set(h1,'color','w','position',[10 100 1900 400]);
ddd=[0 0 0;jet(10)];
orient landscape
RelNatZ = (RelNat - mean(RelNat(:))) ./ std(RelNat(:));
RelGratZ = (RelGrat - mean(RelNat(:))) ./ std(RelNat(:));
% subplot(1,5,1); imagesc( avgProj); colormap gray; axis square; axis off; title('AVG Projection'); freezeColors;
subplot(2,3,1); imagesc( vizMaxProjection(natActMap) ); colormap(ddd); axis square; axis off; title('NS Activity');%caxis([0,700]);
subplot(2,3,2); imagesc( RelNatZ ); colormap(ddd); axis square; axis off; title('NS Rel');%caxis([0, 0.3]);
subplot(2,3,3); imagesc( normRelNat ); colormap(ddd); axis square; axis off; title('NS Rel norm');%caxis([0, 0.3]);
subplot(2,3,4); imagesc( vizMaxProjection(gratActMap) ); colormap(ddd); axis square; axis off; title('Grat Activity');%caxis([0,700]);
subplot(2,3,5); imagesc( RelNGratZ ); colormap(ddd); axis square; axis off; title('Grat Rel');%caxis([0, 0.3]);
subplot(2,3,3); imagesc( normRelGrat ); colormap(ddd); axis square; axis off; title('Grat Rel norm');%caxis([0, 0.3]);

% Create second figure for just the PS plot, passing the axesMod values only.
h2 = figure;
sublot(1,2,1)
plotPowerSpectra(RelGrat, RelNat, axisMod, h2,fname)
title('PS','FontSize',10,'FontWeight','normal');
sublot(1,2,2)
plotPowerSpectra(RelGrat, RelNat, axisMod, h2,fname)
title('PS','FontSize',10,'FontWeight','normal');

%% Save 
if saveFlag
    set (h1, 'PaperOrientation','portrait','PaperPosition',[1 1 14 4]);
    figName1 = [fname(1:end-8), 'RelMap.eps'];
    figName2 = [fname(1:end-8), 'RelMap.jpg'];
    figName3 = [fname(1:end-8), 'PS.eps'];
    figName4 = [fname(1:end-8), 'PS.jpg'];
    
    if saveFlag
        fprintf('Saving figures...\n');
        saveas(h1, figName1, 'eps2c');
        saveas(h1, figName2, 'jpeg');
        saveas(h2, figName3, 'eps2c');
        saveas(h2, figName4, 'jpeg');
    end
    pause(2)
    % close all
end
%% Copy tp relmaps folde for population analysis

if copyFlag
    natFile = [fname(1:end-8), 'RelNat.mat'];
    gratFile = [fname(1:end-8), 'RelGrat.mat'];
    targetDir = 'H:\Dropbox (MIT)\Astrocytes\Data\2pImaging\RelMaps\SingleAstros';
    button = questdlg('Copy files for population analysis?','Copy to RelMaps','Yes');
    
    if strcmp(button, 'Yes')
        [copyStatus, errmsg] = copyfile(figName1,targetDir);
        if copyStatus == 1
            fprintf('%s copied!\n', figName1);
        else
            fprintf('%s NOT copied!\n', figName1);
            disp(errmsg);
        end
        
        [copyStatus, errmsg] = copyfile(figName2,targetDir);
        if copyStatus == 1
            fprintf('%s copied!\n', figName2);
        else
            fprintf('%s NOT copied!....', figName2);
            disp(errmsg);
        end
        
        [copyStatus, errmsg] = copyfile(natFile,targetDir);
        if copyStatus == 1
            fprintf('%s copied!\n', natFile);
        else
            fprintf('%s NOT copied!....', natFile);
            disp(errmsg);
        end
        
        [copyStatus, errmsg] = copyfile(gratFile,targetDir);
        if copyStatus == 1
            fprintf('%s copied!\n', gratFile);
        else
            fprintf('%s NOT copied!....', gratFile);
            disp(errmsg);
        end
%         close all
    end
end
    