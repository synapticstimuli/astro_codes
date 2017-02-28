function plotRelMaps(fname, RelNat, RelGrat)
%% Plots the reliability maps for gratings and ns stims, as well as 
%  the powers spectra analysis. Saves all graphs as eps and jpgs. 
%  From RVR, modified for astros by RG.

%% Process input arguments and input params
if nargin==0
    fprintf('Select data.mat file...\n\n')
    [fname, pathname] = uigetfile('*data.mat','Load data.mat file containing data structure...');
    cd(pathname);
    fprintf('Loading %s...\n\n', fname);
    load (fname, 'data')
    load([pathname fname(1:end-8), 'RelNat']);
    load([pathname fname(1:end-8), 'RelGrat']);
    % Load avg proj image..............................
    try
        projFile = ['AVG_',fname(1:end-9),'.jpg'];
        avgProj = imread(projFile);
    catch
        projFile = uigetfile('*.jpg','Select projection jpg file...');
        avgProj = imread(projFile);
    end
elseif nargin==3
    load (fname, 'data')
    fprintf('Loading data to plot...\n');
end

% Variables 
axisMod         = [2 11];           % Range of x axis for plotPowerSpectra, leave empty for fulle range based on data set
saveFlag        =   1;
plotFlag        =   1;
cropFlags       = [96, 128, 256];   % Sizes indicating files have not been cropped to remove edge effects/
cropSz          =   4;               % Amount to crop from each edge (eg 4 + 1)
cropSt          =   cropSz+1;          
copyFlag        = 0;

%% Plot Rel Maps

% Crop file if necessary
if ismember(size(RelNat,1), cropFlags)
    fprintf('Cropping movie files...\n');
    % Crop image matrix to remove registration edge effects.
    RelNat = RelNat(cropSt:end-cropSz,cropSt:end-cropSz);
    RelGrat = RelGrat(cropSt:end-cropSz,cropSt:end-cropSz);
    avgProj = avgProj(cropSt:end-cropSz,cropSt:end-cropSz);
    dataDFF = data.dFF(cropSt:end-cropSz,cropSt:end-cropSz,:);
    fprintf('Movie files cropped...\n');
else
    dataDFF = data.dFF;
end

% RelNat.............................................
medRelNat=median(median(RelNat));
RelNat(RelNat<medRelNat) = NaN;
minvNat = min(min(RelNat));
maxvNat = max(max(RelNat));
RelNat(isnan(RelNat)) = minvNat-((maxvNat-minvNat)/5);

% RelGrat.............................................
meanRelGrat=mean(mean(RelGrat));
RelGrat(RelGrat<meanRelGrat) = NaN;
minvGrat = min(min(RelGrat));
maxvGrat = max(max(RelGrat));
RelGrat(isnan(RelGrat)) = minvGrat-((maxvGrat-minvGrat)/5);

RelNatZ = (RelNat - mean(RelNat(:))) ./ std(RelNat(:));
RelGratZ = (RelGrat - mean(RelNat(:))) ./ std(RelNat(:));

% Plot first figure containing all plots...............
h1 = figure; 
set(h1,'color','w','position',[10 100 1900 400]);
ddd=[0 0 0;jet(10)];
orient landscape

subplot(1,5,1); 
    imagesc( avgProj); 
    colormap gray; 
    axis square; 
    axis off; 
    title('AVG Projection'); 
    freezeColors;
subplot(1,5,2); 
    imagesc( vizMaxProjection(data.dFF) ); 
    colormap(ddd); 
    axis square; 
    axis off; 
    title('Activity');%caxis([0,700]);
subplot(1,5,3); 
    imagesc( RelNatZ ); 
    colormap(ddd); 
    axis square; 
    axis off; 
    title('NS Map');%caxis([0, 0.3]);
subplot(1,5,4); 
    imagesc( RelGratZ ); 
    colormap(ddd); 
    axis square; 
    axis off; 
    title('Gratings Map');%caxis([0, 0.3]);
subplot(1,5,5);
    plotPowerSpectra(RelGrat,RelNat,axisMod,h1,fname);
    title('PS'); 
    legend('off')

% Plot second figure for just the PS plot, passing the axisMod values only.
h2 = plotPowerSpectra(RelGrat,RelNat,axisMod);
    title(num2str(fname(1:end-18)),'FontSize',10,'FontWeight','normal');

%% Save figures as jpg and esp files

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

%% Copy to relmaps folde for population analysis

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
        close all
    end
end
    