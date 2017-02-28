function [] = RelMapsPop(folder)
%% This will load all the RelGrat and RelNat files within a selected

% directory and extract individual cell/file data for population processing
% Outputs GratRelData, NatRelData, and Figures.

setdlg = 1;
if setdlg
    settings = settingsdlg(...
        'Description'                                       ,'Set parameters:'                      ,...
        'title'                                             ,'Rel Map Population analysis options'  ,...
        {'Select Protocol','protocol'}                      ,{'-NS-','-NSK', '-NS2'}                ,...
        {'Experiment string','expTitle'}                    ,{'','(Pre-Dox)','(Dox-2WK)','(Dox-4WK)',...
                                                             '(Het Pre-IGF1)','(Het Post-IGF1)'     ,...
                                                             '(WT Pre-IGF1)','(WT Post-IGF1)','(Cef7d)'}      ,...
        {'Set YLims','ylims'}                               ,'1 1000'                               ,...
        'separator'                                         ,'Select data to process:'              ,...
        {'Process Rel Map Files','proRelFiles'}             ,[true]                                 ,...
        'separator'                                         ,'Select plots:'                        ,...
        {'All PS and population means','plotAllPS'}         ,[true]                                 ,...
        {'Slopes','plotAllSlopes'}                          ,[true]                                 ,...
        {'Coefficients (box)','plotSlopesBox'}              ,[true]                                 ,...
        {'Coefficients (histogram)','plotSlopesHist'}       ,[false]                                ,...
        {'Area Under the Curve (box)','plotAucBox'}         ,[true]                                 ,...
        {'Area Under the Curve (histogram)','plotAucHist'}  ,[false]                                ,...
        {'Values Table','plotTable'}                        ,[true]                                 ,...
        {'All Maps', 'plotMaps'}                            ,[true]                                 ,...
        'separator'                                         ,'Normalize using blank activity:'      ,...
        {'Include blank activity data?','proActivity'}      ,[false]                                ,...
        {'Load existing activity data?', 'loadActivityData'},[false]                                ,...
        {'Blank Activity Slopes','plotActivitySlopes'}      ,[false]                                ,...
        {'Normalize and plot','plotNormSlopes'}             ,[false]                                ,...
        'separator'                                         ,'Save adn Finish:'                                ,...
        {'Save all variables:','saveVars'}                  ,[true]                                 ,...
        {'Save all figures:','saveFigs'}                    ,[true]                                 ,...s
        {'Close all figures', 'closeAll'}                   ,[false]                                ...
    );
    settings.subPop = 0;
    settings.ylims = str2num(settings.ylims);
else    
    % Files and data to process...............................................
    settings.protocol = '-NS-';
    settings.subPop = 0;
    settings.proRelFiles = 1;        % 1 to process all the rel map files or 0 load processed values
    settings.proActivity = 0;        % 1 to include activity maps or 0 skip this analysis
    settings.loadActivityData = 0;   % 1 to load pre processed activity data or 0 process files
    
    % Choose what to save......................................................
    settings.saveVars = 1;
    settings.saveFigs = 1;
    
    % Choose which figuresd to generate........................................
    settings.plotAllPS = 1;          % Plot figure with all PS and populatuon means subplots
    settings.plotAucHist = 0;        % Plot AuC Histogram
    settings.plotAucBox = 1;         % Plot Auc boxScatter plot
    settings.plotAllSlopes = 0;      % Plot all slopes(gradients)
    settings.plotActivitySlopes = 0; % Plot slopes of blank activity frames
    settings.plotSlopesBox = 1;      % Plot bo scatter ploat of slopes
    settings.plotSlopesHist = 0;     % Plot histogram of slopes
    settings.plotNormSlopes = 0;     % Plot normalized values slopes
    settings.plotTable = 1;          % Plot table with all coefficients
    settings.plotMaps = 1;
    % Other input params
    settings.expTitle = '(Het Pre-IGF1)';     % Figure window name modifier, eg '(Post-Dox 4w(K2))'
    settings.ylims = [1 1000];
end

%% Select directory
if nargin==0
    folder = uigetdir;
    filelist = dir([folder filesep '*' settings.protocol '*']);
    nfiles = length(filelist);
    if ~isempty(nfiles)
        fprintf('%1.0f NS and Grat files found...\n',nfiles);
        pause(1)
    elseif isempty(nfiles)
        fprintf('%1.0f files found...\n',nfiles);
        pause(5)
    end
        
end

dataFolder = [folder filesep 'data' filesep];
figFolder = [folder filesep 'figures' filesep];
%% Load reliability files
if settings.proRelFiles
    close all
    
    fprintf('Processing reliability files...\n');
    
    cropFlags = [96, 128, 256]; %Set sizes to check for uncropped files that may have edge effects.
    cropSz = 4;
    
    gratIndx = 0;
    natIndx = 0;
    maxSz = 65;
    natRelData = struct('fname',{},'data',{});
    gratRelData = struct('fname',{},'data',{});
    
    for i=3:length(filelist)
        
        %load RelNat files
        if ~isempty (strfind(filelist(i).name,'RelNat'))
            disp(filelist(i).name);
            natIndx = natIndx + 1;
            load ([folder, '/',filelist(i).name])
            
            % Crop rel maps to remove edge effects
            if ismember(size(RelNat,1), cropFlags)
                RelNat = RelNat(cropSz:end-cropSz, cropSz:end-cropSz);
                fprintf('%s cropped...\n',filelist(i).name);
            end
            %         pause(1)
            
            natRelData{natIndx}.fname = filelist(i).name;
            natRelData{natIndx}.data = RelNat;
            AX = fftshift(fft2(RelNat));
            ASNat{natIndx} = abs(AX).^2;
            PSNat{natIndx} = rotavg(ASNat{natIndx});
            PSSize = length(PSNat{natIndx});
            %         if maxSz < PSSize
            %             maxSz = PSSize;
            %         end
            if PSSize < maxSz
                PSNat{natIndx} = [PSNat{natIndx}; NaN(maxSz-length(PSNat{natIndx}),1)];
            end
            clear RelNat
            
            %load GratNat files
        elseif ~isempty (strfind(filelist(i).name,'RelGrat'))
            disp(filelist(i).name);
            gratIndx = gratIndx + 1;
            load ([folder, '/',filelist(i).name])
            
            % Crop rel maps to remove edge effects
            if ismember(size(RelGrat,1), cropFlags)
                RelGrat = RelGrat(cropSz:end-cropSz, cropSz:end-cropSz);
                fprintf('%s cropped...\n',filelist(i).name);
            end
            %         pause(1)
            gratRelData{gratIndx}.fname = filelist(i).name;
            gratRelData{gratIndx}.data = RelGrat;
            AX = fftshift(fft2(RelGrat));
            ASGrat{gratIndx} = abs(AX).^2;
            PSGrat{gratIndx} = rotavg(ASGrat{gratIndx});
            PSSize = length(PSGrat{gratIndx});
            %         if maxSz < PSSize
            %             maxSz = PSSize;
            %         end
            if PSSize < maxSz
                PSGrat{gratIndx} = [PSGrat{gratIndx}; NaN(maxSz-length(PSGrat{gratIndx}),1)];
            end
            clear RelGrat
        end
    end
    
    fprintf ('%2.0f grat files found \n',length(gratRelData));
    fprintf ('%2.0f nat files found \n',length(natRelData));
    
    if settings.saveVars
        save ([folder filesep 'data' filesep 'GratRelData.mat'],'gratRelData', 'PSGrat', 'ASGrat');
        save ([folder filesep 'data' filesep 'NatRelData.mat'],'natRelData','PSNat','ASNat');
    end
else
    fprintf('Loading reliability data files...\n');
    load([dataFolder 'GratRelData.mat']);
    load([dataFolder 'NatRelData.mat']);
    load([dataFolder 'PSdata.mat']);
end
%% Load activity files files

if settings.proActivity
    if settings.loadActivityData
        fprintf('Loading activity data files...\n');
        load([dataFolder 'GratActData.mat'])
        load([dataFolder 'NatActData.mat'])
    else
        fprintf('Loading activity files...\n');
        natAct = dir([folder filesep '*NatImages.mat']);
        numnatAct = length(natAct);
        fprintf('%1.0f NS activity files found.\n',numnatAct);
        gratAct = dir([folder filesep '*Gratings.mat']);
        numgratAct = length(gratAct);
        fprintf('%1.0f Grat activity files found.\n',numgratAct);
        
        for i = 1: numgratAct
            % Extract mean and max activity projections for every experiment file
            natActData{i}.fname = natAct(i).name;
            load([folder filesep natActData{i}.fname]);
            currNatdata = NatImages;
            for j = 1:6
                curr = currNatdata{j};
                curMax(:,:,j) = max(curr,[],3);
            end
            natActData{i}.max = max(curMax,[],3);
            natActData{i}.mean = mean(curMax,3);
            clearvars NatImages currNatdata curMax curr
            
            % Activity of grating responses
            gratActData{i}.fname = gratAct(i).name;
            load([folder filesep gratActData{i}.fname]);
            currGratdata = Gratings;
            for j = 1:6
                curr = currGratdata{j};
                curMax(:,:,j) = max(curr,[],3);
            end
            gratActData{i}.max = max(curMax,[],3);
            gratActData{i}.mean = mean(curMax,3);
            clearvars Gratings currGratdata curMax curr
            fprintf('%1.0f of %1.0f activity files processed...\n',i, numnatAct);
        end
        
        if settings.saveVars
            fprintf('Saving activity files...\n');
            save ([folder filesep 'data' filesep 'GratActData.mat'],'gratActData');
            save ([folder filesep 'data' filesep 'NatActData.mat'],'natActData');
        end
    end
end
%% Plot population powerspectra

% Load population PS data and calculate means and AUC.
xNat = cell2mat(PSNat);
xNat = xNat(3:12,:);
xNatMean = nanmean(xNat,2);
xNatSEM  = nanstd(xNat,[],2)/sqrt(length(xNat));
aocNat = trapz(xNat);

xGrat = cell2mat(PSGrat);
xGrat = xGrat(3:12,:);
xGratMean = nanmean(xGrat,2);
xGratSEM  = nanstd(xGrat,[],2)/sqrt(length(xGrat));
aocGrat = trapz(xGrat);

if settings.plotAllPS    % Plot all traces and population means
    f1 = figure; hold on; 
    subplot(1,2,1)
    hold on; 
    loglog(xNat,'Color',[.6 .6 .6],'LineWidth',0.5);
    loglog(xGrat,'Color',[1 .4 .4],'LineWidth',0.5);
    % loglog(xNatMean,'k','LineWidth',6);
    % loglog(xGratMean,'r','LineWidth',6);
    title('PS - All Cells');
    axis square; box off;
    set(gca,'xscale','log','yscale','log','YLim',settings.ylims);
%     settings.ylims = get(gca,'Ylim');
    
    subplot(1,2,2)
    hold on; 
    loglog(xNatMean,'k','LineWidth',6);
    loglog(xGratMean,'r','LineWidth',6);
    loglog(xNatMean + xNatSEM,'k','LineWidth',1);
    loglog(xNatMean - xNatSEM,'k','LineWidth',1);
    loglog(xGratMean + xGratSEM,'r','LineWidth',1);
    loglog(xGratMean - xGratSEM,'r','LineWidth',1);
    title('Population PS Mean');
    axis square; box off;
    
    set(gca,'xscale','log','yscale','log','YLim', settings.ylims);
    set(f1,...
        'color'         ,'w'                          ,...
        'Name'          ,['All_PS' settings.expTitle]           ,...
        'Position'      ,[100 100 800 300],...
        'NumberTitle'   ,'off'    );
end
    
if settings.plotAucHist  % Plot histogram for area under the curve
    f2 = figure; hold on;
    h3 = histogram(aocGrat,10);
    h4 = histogram(aocNat,10);
    h3.Normalization = 'probability';
    h3.BinWidth = 25;
    h4.Normalization = 'probability';
    h4.BinWidth = 25;
    legend('Grat','Nat');
    title('Area Under the Curve');
    
    set(f2,...
        'color'         ,'w'                    ,...
        'Name'          ,['AuC_Hist' settings.expTitle]  ,...
        'NumberTitle'   ,'off'                  );
end

if settings.plotAucBox   % Plot box + scatter plot of area under the curve
    f3 = figure;
    boxScatter([aocNat' aocGrat'], f3);
    title ('Area Under the Curve');
    
    set(f3,...
        'color'         ,'w'                            ,...
        'Name'          ,['AuC_Box' settings.expTitle]          ,...
        'NumberTitle'   ,'off'                          );
end
%% Plot population PS slopes

if settings.plotAllSlopes
    gratFiles = length(gratRelData);
    natFiles = length(natRelData);
    logfreqNat = cell(1,natFiles);
    logPfNat = cell(1,natFiles);
    logfreqGrat = cell(1,gratFiles);
    logPfGrat = cell(1,gratFiles);
    
    %   Plots the coefficients/gradients for the NS of each cell on one subplot
    f4 = figure; hold on
    subplot(1,2,1)
    for j = 1:natFiles
        xN = natRelData{j}.data;
        [logfreqNat{j}, logPfNat{j}, natSlopes(j)] = plotPowerSpectrumRVR(xN,f4);
        clc
    end
    ax = gca;
    ax.Title.String = 'Nat Scenes';
    ax.Title.FontWeight = 'normal';
    settings.ylims = ax.YLim;
    hold off;
    
    %   Plots the coefficients/gradients for the Gratings of each cell on one subplot
    subplot(1,2,2)
    for j = 1:gratFiles
        xG=gratRelData{j}.data;
        [logfreqGrat{j}, logPfGrat{j}, gratSlopes(j)] = plotPowerSpectrumRVR(xG,f4);
        clc
    end
    ax = gca;
    ax.Title.String = 'Gratings';
    ax.Title.FontWeight = 'normal';
    ax.YLim = settings.ylims;
    hold off;
    set(f4,...
        'color'         ,'w'                            ,...
        'Name'          ,['All_Slopes' settings.expTitle]       ,...
        'Position'      ,[200 200 1000 300]             ,...
        'NumberTitle'   ,'off'                          );
    
    if settings.saveVars
        save ([folder filesep 'data' filesep 'GratRelData.mat'],'gratSlopes', 'aocNat', '-append');
        save ([folder filesep 'data' filesep 'NatRelData.mat'],'natSlopes','aocGrat', '-append');
    end
else
    if settings.plotActivitySlopes
        %   Plots the coefficients/gradients for the blank stim frames of 
        %   the NS of each cell on one subplot
        f5 = figure; hold on
        subplot(1,2,1)
        for j = 1:natFiles
            xNA = natActData{j}.max;
            [logfreqNatAct{j}, logPfNatAct{j}, natActSlopes(j)] = plotPowerSpectrumRVR(xNA,f5);
            clc
        end
        ax = gca;
        ax.Title.String = 'Nat Scenes Normalized';
        ax.Title.FontWeight = 'normal';
        settings.ylims = ax.YLim;
        hold off;
        
        subplot(1,2,2)
        %   Plots the coefficients/gradients for the blank stim frames of 
        %   the Gratings of each cell on one subplot
        for j = 1:gratFiles
            xGA = gratActData{j}.max;
            [logfreqGratAct{j}, logPfGratAct{j}, gratActSlopes(j)] = plotPowerSpectrumRVR(xGA,f5);
            clc
        end
        ax = gca;
        ax.Title.String = 'Gratings Normalized';
        ax.Title.FontWeight = 'normal';
        ax.YLim = settings.ylims;
        hold off;
        set(f5,...
            'color'         ,'w'                                ,...
            'Name'          ,['All_Slopes_Norm' settings.expTitle]      ,...
            'Position'      ,[200 200 1000 300]                 ,...
            'NumberTitle'   ,'off'                              );
        
        % Generate normalized slopes (by activity coeeficients)
        natSlopesNorm = natSlopes ./ natActSlopes;
        gratSlopesNorm = gratSlopes ./gratActSlopes;
        
        if settings.saveVars
            save ([folder filesep 'data' filesep 'GratActData.mat'],'gratActSlopes', '-append');
            save ([folder filesep 'data' filesep 'NatActData.mat'],'natActSlopes', '-append');
        end
        
    end
end 
if settings.plotSlopesBox    % Plot box + scatter plot of slopes
    f6 = figure; hold on;
    boxScatter([abs(natSlopes') abs(gratSlopes')], f6,{'NS' 'Gratings'});
    set(gca, 'YLim', [0 4]);
    title ('PS Slopes');
    set(f6,...
        'color'         ,'w'                            ,...
        'Name'          ,['Slopes_Box' settings.expTitle]       ,...
        'NumberTitle'   ,'off'                          );
end
if settings.plotNormSlopes   % Plot box + scatter plot of normalized slopes
    f7 = figure; hold on;
    boxScatter([abs(natSlopesNorm') abs(gratSlopesNorm')], f7,{'NS' 'Gratings'});
    title('PS Slopes Normalized');
    set(f7,...
        'color'         ,'w'                                ,...
        'Name'          ,['Slopes_Norm_Box' settings.expTitle]      ,...
        'NumberTitle'   ,'off'                              );
end
if settings.plotSlopesHist   % Plot histogram of slopes
    f8 = figure; hold on;
    h9a = histogram(natSlopes,20);
    h9b = histogram(gratSlopes,20);
    legend('Grat','Nat');
    set(f8,...
        'color'         ,'w'                                ,...
        'Name'          ,['Slopes_Hist' settings.expTitle]      ,...
        'NumberTitle'   ,'off'                              );
    hold off;
end
if settings.subPop           % Plot sub- and supra-mean populations
    f9 = subPopulationAnalysis(natSlopes, gratSlopes,natRelData,gratRelData);
    set(f9,...
        'color'         ,'w'                                ,...
        'Name'          ,['SubPop Analysis' settings.expTitle]      ,...
        'NumberTitle'   ,'off'                              );
end
%% Create table with slope coefficients
if settings.plotTable
    f10 = figure;
    for i = 1:length(gratRelData)
        tableData{i,1} = gratRelData{i}.fname(1:end-12);
        tableData{i,2} = abs(natSlopes(i));
        tableData{i,3} = abs(gratSlopes(i));
        tableData{i,4} = abs(aocNat(i));
        tableData{i,5} = abs(aocGrat(i));
        if settings.plotNormSlopes
            tableData{i,6} = abs(natSlopesNorm(i));
            tableData{i,7} = abs(gratSlopesNorm(i));
        end
    end
    tableData = sortcell(tableData,2);
    
    colNames = {'File Names' 'NS Coef' 'Grat Coef' 'Nat AUC' 'Grat AUC' ' NS nCoef' 'Grat nCoef'};
    colForm = { 'char', 'numeric', 'numeric','numeric','numeric','numeric','numeric'};
    
    t = uitable(f10, 'Data', tableData, 'ColumnName', colNames, 'ColumnFormat', colForm);
    
    % Set width and height
    t.Position(3) = 950;
    t.Position(4) = 750;
    t.ColumnWidth={350,'auto','auto','auto','auto','auto','auto'};
    
    set(f10,...
        'color'         ,'w'                                        ,...
        'Name'          ,['All Values' settings.expTitle]                   ,...
        'Position'      ,[100 100 t.Position(3)+50 t.Position(4)+50],...
        'NumberTitle'   ,'off'                                      );
    
    gradNatMean = mean(abs(natSlopes));
    gradGratMean = mean(abs(gratSlopes));
    
    fprintf('NS mean: %2.2f\n', gradNatMean);
    fprintf('Grat mean: %2.2f\n', gradGratMean);
    
    
    
    if settings.saveVars
        colnames = get(t,'columnname')';     % t is the handle of your table
        data=get(t,'data');
        data=[colnames(1:5);data];
        fid = fopen([folder filesep 'figures' filesep 'AllValues' settings.expTitle '.csv'],'w');
        fprintf(fid,'%s, %s, %s, %s, %s\n',data{1,1:5});
        %         fprintf(fid,'%s\n',data{2:end,1})
        for i = 2:length(data)
            fprintf(fid,'%s, %f, %f, %f, %f\n',data{i,1:5});
        end
    end
end
%% Plot all maps
if settings.plotMaps
    mapFiles = dir('TSeries*.jpg');
    numImg = length(mapFiles);
    figs = ceil(numImg/6);
    ct = 1;
    for i = 1:figs
        f11(i) = figure; hold on;
        for j = 1:6
            if ct < numImg +1
                subplot_tight(3,2,j,[0.05 0.001]);
                curName = mapFiles(ct).name;
                curImg = imread(curName);
                imagesc(curImg);
                title(curName(1:end-11),'FontSize', 8, 'FontWeight', 'normal');
                axis off;
                set(gcf,...
                    'color'         ,'w'                            ,...
                    'Name'          ,['All Maps' settings.expTitle '-0' num2str(i)]   ,...
                    'Position'      ,[50 50 1200 700]               ,...
                    'NumberTitle'   ,'off'                          );
                ct = ct + 1;
            end
        end
    end
end
%% Save variables
if settings.saveVars
    save([folder filesep 'data' filesep 'PSdata.mat'],'logfreqNat','logPfNat','natSlopes','logfreqGrat','logPfGrat','gratSlopes','aocGrat','aocNat');
end

%% Save figures
if settings.saveFigs
    fprintf('Saving figures...\n');
    t = get(0,'children');
    cd(['..' filesep])
    cd('figures');
    for i=1:length(t)
        saveas(t(i),[t(i).Name '.jpg'],'jpeg')
        saveas(t(i),[t(i).Name '.eps'],'eps2c')
    end
end
cd(folder);

if settings.closeAll
    close all
end

%% 
