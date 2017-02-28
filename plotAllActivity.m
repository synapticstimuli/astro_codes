function plotAllActivity(fname)

% plots two figures given the data file for an NS experiment. The first
% contains 6 subplots: avg activity, NS map, grat map, dff max proj, ns act
% max proj, and grat activity max proj. 
% The second plots the individual trial max dff activity for ns and grat.

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
elseif nargin == 1
    fprintf('Loading files...\n\n');
    load([fname(1:end-8), 'RelNat']);
    load([fname(1:end-8), 'RelGrat']);
    load([fname(1:end-8), 'NatImages']);
    load([fname(1:end-8), 'Gratings']);
elseif nargin==3
    load (fname, 'data')
    fprintf('Loading data to plot...\n');
end

%% Input params

cropFlags       = [96, 128, 256];   % Sizes indicating files have not been cropped to remove edge effects/
cropSz          =   7;               % Amount to crop from each edge (eg 4 + 1)
cropSt          =   cropSz+1;  
%% Crop files if necessary

fprintf('Cropping movie files...\n');
% Crop image matrix to remove registration edge effects.
RelNat = RelNat(cropSt:end-cropSz,cropSt:end-cropSz);
RelGrat = RelGrat(cropSt:end-cropSz,cropSt:end-cropSz);
for i =1:6;
    NatImages{i} = NatImages{i}(cropSt:end-cropSz,cropSt:end-cropSz,:);
    Gratings{i} = Gratings{i}(cropSt:end-cropSz,cropSt:end-cropSz,:);
end
%     avgProj = avgProj(cropSt:end-cropSz,cropSt:end-cropSz);
dataDFF = data.dFF(cropSt:end-cropSz,cropSt:end-cropSz,:);
fprintf('Movie files cropped...\n');

%% Get response activity
natActMap=[];
for i=1:length(NatImages)
    natActMap = cat(3,natActMap, NatImages{i});
end
gratNatMap=[];
for i=1:length(Gratings)
    gratNatMap = cat(3,gratNatMap, Gratings{i});
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

%Load avg proj image
    try
        projFile = ['AVG_',fname(1:end-9),'.jpg'];
        avgProj = imread(projFile);
    catch
        projFile = uigetfile('*.jpg','Select projection jpg file...');
        avgProj = imread(projFile);
    end

if nargin == 0
    f1 = figure; 
    set(gcf,'color','w','position',[10 100 1800 800]);
    ddd=[0 0 0;jet(10)];
    orient landscape
    
    subplot(231); imagesc( avgProj); colormap gray; axis square; axis off; title('AVG Projection'); freezeColors;
    subplot(2,3,2); imagesc( RelNat ); colormap(ddd); axis square; axis off; title('NS Rel');
    subplot(2,3,3); imagesc( RelGrat ); colormap(ddd); axis square; axis off; title('Grat Rel');
    subplot(234); imagesc( vizMaxProjection(dataDFF) );colormap(ddd); axis square; axis off; title('Activity');
    subplot(235); imagesc( vizMaxProjection(natActMap) ); colormap(ddd); axis square; axis off; title('NS Activity');
    subplot(236); imagesc( vizMaxProjection(gratNatMap) ); colormap(ddd); axis square; axis off; title('Grat Activity');
    shg
    
    set(f1,...
        'color'         ,'w'                            ,...
        'Name'          ,['All_Maps-' fname]           ,...
        'Position'      ,[200 200 1800 800]             ,...
        'NumberTitle'   ,'off'                          );
end

f2 = figure; 
set(gcf,'color','w','position',[10 100 1800 800]);
ddd=[0 0 0;jet(10)];
orient landscape
for j = 1:6
    currnat = NatImages{j};
    subplot_tight(2,6,j,[0.001 0.01]);
    imagesc(vizMaxProjection(currnat));colormap(ddd); axis square; axis off; title(['NS',num2str(j)]);
    currgrat = Gratings{j};
    subplot_tight(2,6,j+6,[0.001 0.01]);
    imagesc(vizMaxProjection(currgrat));colormap(ddd); axis square; axis off; title(['Grat',num2str(j)]);
end
set(f2,...
        'color'         ,'w'                            ,...
        'Name'          ,['NS_Grat_Maps-' fname]           ,...
        'Position'      ,[200 200 1800 800]             ,...
        'NumberTitle'   ,'off'                          );

%% Save figures
% fprintf('Saving figures...\n');
% saveas(f1, [fname(1:end-8) 'ActivityMaps-01'],'jpg');
% saveas(f1, [fname(1:end-8) 'ActivityMaps-01'],'eps2c');
% saveas(f2, [fname(1:end-8) 'ActivityMaps-02'],'jpg');
% saveas(f2, [fname(1:end-8) 'ActivityMaps-02'],'eps2c');
%     
 fprintf('Saving figures...\n');
    t = get(0,'children');
    for i=1:length(t)
        saveas(t(i),[t(i).Name '.jpg'],'jpeg')
        saveas(t(i),[t(i).Name '.eps'],'eps2c')
    end
    