function roiPlot_exp(fname)
% Select rois.mat file to run...
if nargin ==0
    clearvars 
    close all
    clc
    roisFile = uigetfile('*-rois.mat','Select *rois.mat file... ');
    try
        load(roisFile);
    catch
        fprintf('*** No file loaded! Script aborted... ***\n');
        return
    end
elseif nargin ==1
    load(fname);
end

%% Begin ui settings input


settings.plotAllTraces  = 0;
settings.plotNMTraces   = 0;
settings.plotGratTraces = 0;
settings.plotIndOri     = 0;
settings.plotSummaryTC  = 0;
settings.shiftFlag      = 1;
settings.subROIs        = 'All';
settings.filtRoiFlag    = 0;
settings.plotOverlay    = 1;
settings.roiType        = 'Rectangular';
settings.plotOSI        = 1;
settings.plotPO         = 0;
settings.saveFigs       = 1;
settings.closeAll       = 0;

% [settings, button] = settingsdlg(...
%     'Description'                                       ,'Set parameters:'                      ,...
%     'title'                                             ,'ROI NM and Ori analysis options'      ,...
%     'separator'                                         ,'Select plots:'                        ,...
%     {'Plot all traces?','plotAllTraces'}                ,[false]                                ,...
%     {'Plot NM traces and raster plots?','plotNMTraces'} ,[false]                                ,...
%     {'Plot gratings and raster plots?','plotGratTraces'},[false]                                ,...
%     {'Plot figure with individual oris?','plotIndOri'}  ,[false]                                ,...
%     {'Plot summary TCs?','plotSummaryTC'}               ,[true]                                 ,...
%     {'Shift TC axes to center at 0°?','shiftFlag' }     ,[true]                                 ,...
%     'separator'                                         ,'Select ROIs:'                         ,...
%     {'Select ROIs to run','subROIs'}                    ,'All'                                  ,...
%     {'Plot only filtered ROIs?','filtRoiFlag'}          ,[true]                                 ,...
%     'separator'                                         ,'Select overlay plots:'                ,... 
%     {'Plot overlay of rois?','plotOverlay'}             ,[true false]                            ,...
%     {'ROI overlay type?','roiType' }                    ,{'Circular','Rectangular'}             ,...
%     {'Plot OSI Overlay?','plotOSI'}                      ,[true]                                 ,...
%     {'Plot PO Overlay?','plotPO'}                        ,[true]                                 ,...
%     'separator'                                         ,'Save and Finish:'                     ,...
%     {'Save all figures?','saveFigs'}                    ,[true]                                 ,...
%     {'Close all figures?', 'closeAll'}                  ,[false]                                ...
%     );
% 
% if strcmp(button,'cancel') || isempty(button)
%     fprintf('\n*****************  No settings were entered. Script stopped.  *********************\n\n\n');
%     return
% end

%% Define variables
% Trial specific structures...
NumSamplsNew    = rois.stimDur .* rois.newFrameRate;                                       % Total number of frames for entire experiment
TimeEachTrial   = (Params.repMovies * (Params.durStim_ns + Params.durBlank_nm))...      % Duration (s) of each trial
                  + (Params.numOris * (Params.durStim_grats + Params.durBlank_grats));
SamplsEachTrial = TimeEachTrial * rois.newFrameRate;                                    % Number of frames for each trial
NumTrials       = NumSamplsNew ./ SamplsEachTrial;                                      % Calculation of num trials (see also Params.repTrials)

% Stimulus-specific trial structures... 
durSec_NM       = Params.durStim_ns + Params.durBlank_nm;                                   % Duration (s) of NS (ON + OFF)
durFrames_NM    = durSec_NM * rois.newFrameRate;                                       % Num frames of NS stim (ON + OFF)
durTrials_NM    = durFrames_NM * Params.repMovies;                                     % Total length of NS (ON + OFF) per trial
numTrials_NM = Params.repMovies * NumTrials;

durSec_Grats    = Params.durStim_grats + Params.durBlank_grats;                        % Duration (s) of Grats (ON + OFF)
durFrames_Grats = durSec_Grats * rois.newFrameRate;                                 % Num frames of Grats stim (ON + OFF)
durTrials_Grats = durFrames_Grats * Params.numOris;
numTrials_Grats = Params.numOris * Params.repTrials;                                % Total numer of gratings repetitions

msgCt = zeros(1,6);
msgCt(1) =  1;

%% Modify roi labels and indices for pre-selected rois
if  ~strcmp(settings.subROIs,'All')
    curROIs = str2num(settings.subROIs);
    NumROIs     = length(curROIs);
    rois.Label = rois.Label(curROIs);
    rois.filtROI = curROIs;
%     settings.filtRoiFlag = 0;
else
    NumROIs = length(rois.Label);
    curROIs = rois.Label;
end

%% Plot all and average traces by stimulus and trials
figCt = 1;

if settings.plotAllTraces
    maxSubPlots = 5;
    subPlotIndx = 0;
    totalCt = NumROIs;
    
    % Create time-based x-axes
    numframes_NM = durSec_NM * Params.repMovies;
    numFrames_Grats = durSec_Grats * Params.numOris;
    xAx_NM = linspace(0, numframes_NM,durTrials_NM);                           % For all NM frames
    xAx_Grats = linspace(0, numFrames_Grats, durTrials_Grats);
   
    for i = 1:NumROIs
        f(figCt) = figure(figCt);
        
        % Plot NS traces
        subPlotIndx = subPlotIndx +1;
        subplot (maxSubPlots, 2, subPlotIndx);
        hold on;
        curROI = reshape(rois.NatMovies(curROIs(i),:), durTrials_NM,[])';
        meanTrace = mean(curROI);
        plot(xAx_NM, curROI','k');
        plot(xAx_NM, meanTrace,'r', 'LineWidth',3);
        set(gca,'XLim', [0 numframes_NM]);
        title(['NS ROI ', num2str(rois.Label(i))]);
        
        % Plot Grat traces
        subPlotIndx = subPlotIndx +1;
        subplot (maxSubPlots, 2, subPlotIndx);
        hold on;
        curROI = reshape(rois.Gratings(curROIs(i),:), durTrials_Grats, [])';
        meanTrace = mean(curROI);
        plot(xAx_Grats, curROI','k');
        plot(xAx_Grats, meanTrace,'r', 'LineWidth',3);
        set(gca,'XLim', [0 numFrames_Grats]);
        title(['Grat ROI ', num2str(rois.Label(i))]);
              
        % Reset subplots and figure counts
        totalCt = totalCt-1;
        if totalCt ~= 0
            if rem(i,maxSubPlots) == 0
                set(gcf                                                 ,... 
                'Position',[(0+(10*figCt)) (0+(10*figCt)) 1300 800] ,...
                'Color', 'w');
                suptitle([rois.fname(1:end-9),' (pg.',num2str(figCt),')']);
                
                subPlotIndx = 0;
                figCt = figCt + 1;
            end
        elseif totalCt == 0
            set(gcf                                                 ,... 
                'Position',[(0+(10*figCt)) (0+(10*figCt)) 1300 800] ,...
                'Color', 'w');
                suptitle([rois.fname(1:end-9),' (pg.',num2str(figCt),')']);
        end
        
    end
    
%     % Save figures...
%     if settings.saveFigs
%         for i = 1:length(f)
%             (f(i));
%             saveName = [rois.fname(1:end-9) '-allTraces-', num2str(i)];
%             saveas(gcf,saveName,'eps2c');
%             saveas(gcf,saveName,'jpeg');
%             fprintf('All traces figure saved!\n');
%             msgCt(2) =  1;
%         end
%      end
end

%% Plot NM traces

if settings.plotNMTraces    
    % Set stimulus box params
    stimON = Params.durBlank_nm;
    stimEnd = durSec_NM;
    
    % Set subplot starting indices
    maxSubPlots = 5;
    subPlotIndx = 0;
    
    % Check existing figures and create figure
    if ishandle(figCt)
        figCt = figCt + 1;
    end
    pgCt = 1;
    
    for ROI = 1:NumROIs;
    
        ff(figCt) = figure(figCt); 
        hold on; 
        set(gcf,'color','w');
    
        resp_NM = squeeze(rois.NatMovies(curROIs(ROI),:,:) );
        resp_NM = resp_NM(:);
        respTr_NM = reshape( resp_NM', [durFrames_NM, Params.repTrials.*Params.repMovies] );
        maxPeak = mean(findpeaks(resp_NM,'MinPeakHeight', 50,'MinPeakDistance',100));
        respTr_NM_norm = respTr_NM./max(respTr_NM(:));
        
        % Create time-based x-axes
        nsTime = length(resp_NM)/rois.newFrameRate;
        xAx_NM = linspace(0, nsTime, length(resp_NM));                      % For all NM frames
        xvals = linspace(0,durSec_NM,durFrames_NM);
    
        % Plot NM trace for al movie reps
        subPlotIndx = subPlotIndx +1;
        subplot_tight(maxSubPlots, 4, [subPlotIndx subPlotIndx+1],[0.1 0.05]);
        hold on;
        plot(xAx_NM, resp_NM, 'linewidth', 1.5, 'color','k'); hold on;
        box off;
        set(gca,'fontname','helvetica','tickdir','out');
        xlim([0,360]);
        ylim([0, max(resp_NM)]);
        for t = 0:15:360
            line([t,t],[0, max(resp_NM)],'linestyle','--','color','r','linewidth',0.5);
            line([t+2,t+2],[0, max(resp_NM)],'linestyle','--','color','b','linewidth',0.5);
        end;
        title(['ROI #',num2str(rois.Label(ROI))]);  
        drawnow;
        
        % Plot NM Raster
        subPlotIndx = subPlotIndx + 2;
        subplot_tight(maxSubPlots, 4, subPlotIndx,[0.1 0.05]);
        
        imagesc(xvals, 1:numTrials_NM, respTr_NM_norm' ); hold on;
        line( [2,2], [1,25],'linewidth', 3, 'color','k');
        colormap( b2r(0,1) ); caxis([0.1,1]);
        axis square; box off; axis xy;
        set(gca,'fontname','helvetica','tickdir','out');
        title(['Rel = ',num2str(rois.Reliability(curROIs(ROI)))]);
        
        % Plot NS ON/OFF Trace
        subPlotIndx = subPlotIndx + 1;
        subplot_tight(maxSubPlots, 4, subPlotIndx,[0.1 0.05]);
        stdshadenew(xvals, mean(respTr_NM,2)', std(respTr_NM,[],2)'./sqrt(numTrials_NM), 0.5, 'k' ); hold on;
        line( [2,2], [1,300],'linewidth', 2, 'color','r');
        axis square; box off;
        set(gca, 'XLim',[0 durSec_NM],...
            'YLim', [0 maxPeak],...
            'fontname','helvetica',...
            'tickdir','out');
        title(['Max = ',num2str(maxPeak)]);
        set(gcf,'Position', [200 50 800 935]);
        suptitle([rois.fname(1:end-9),' (pg.',num2str(pgCt),')']);
        
        % Reset subplots
        if rem(ROI,maxSubPlots) == 0
            subPlotIndx = 0;
            figCt = figCt + 1;
            pgCt = pgCt + 1;
        end
    end
    
    % Save figures...
%     if settings.saveFigs
%         for i = 1:length(ff)
%             (ff(i));
%             saveName = [rois.fname(1:end-9) '-NMplots-', num2str(i)];
%             saveas(gcf,saveName,'eps2c');
%             saveas(gcf,saveName,'jpeg');
%             fprintf('NM traces and raster plots saved!\n');
%             msgCt(3) =  1;
%         end
%     end
end;

%% Plot Ori and TC

if settings.plotIndOri
    % Set subplot starting indices
    maxSubPlots = 5;
    subPlotIndx = 0;
    figCt = figCt + 1;
    maxRows = ceil(NumROIs/5);
    
   
    
    % Check existing figures and create figure
    if ishandle(figCt)
        figCt = figCt + 1;
    end
    pgCt = 1;
    
    Angles = linspace(0,330, 20000);
    AngRad = Angles.*(pi/180);
    
    for ROI = 1:NumROIs
        
        resp_Grat = squeeze(rois.Gratings(curROIs(ROI),:,:) );
        resp_Grat = resp_Grat(:);
        respTr_Grat     = reshape( resp_Grat', [durFrames_Grats, numTrials_Grats] );
        respTr_Grat_norm = respTr_Grat./max(respTr_Grat(:));
        
        % Create time-based x-axes
        gratTime = length(resp_Grat)/rois.newFrameRate;
        xAx_Grats = linspace(0, gratTime, length(resp_Grat));
        xvals = linspace(0,durSec_Grats,durFrames_Grats);
        
        for ori = 1:Params.numOris
            ct = 0;
            for ind = ori : Params.numOris : 72
                ct = ct+1;
                Resp{ROI}{ori}(ct,:) = respTr_Grat( :, ind)';
            end;
            
        end;
        
        % Plot gratings trace for all orientaion reps
        if settings.plotGratTraces
            fff(ROI) = figure(figCt); 
            hold on; 
            set(gcf,'color','w');
            
            subPlotIndx = subPlotIndx +1;
            subplot(maxSubPlots, 1, subPlotIndx);
            hold on;
            plot(xAx_Grats, resp_Grat, 'linewidth', 1.5, 'color','k'); hold on;
            box off;
            set(gca,'fontname','helvetica','tickdir','out');
            xlim([0, 432]); ylim([0, max(resp_Grat)]);
            for t = 0:72:432
                line([t,t],[0, max(resp_Grat)],'linestyle','--','color','r','linewidth',0.5);
                line([t+4,t+4],[0, max(resp_Grat)],'linestyle','--','color','b','linewidth',0.5);
            end;
            ylabel(['ROI #', num2str(rois.Label(ROI))]);
            suptitle([rois.fname(1:end-9),' (pg.',num2str(figCt),')']);
        end
        
        % Reset subplots
        if rem(ROI,maxSubPlots) == 0
            subPlotIndx = 0;
            figCt = figCt + 1;
        end
        
        
        ffff(ROI) = figure(100+ROI); hold on;
        set(gcf,...
            'color','w',...
            'Position', [(1+(ROI*20)) 1 665 981],...
            'MenuBar', 'none'...
            );
        
        for ori = 1:12
            subplot_tight(9,4,ori,[0.03 0.03]);
            imagesc(xvals, 1:6, Resp{ROI}{ori}./max( Resp{ROI}{ori}(:)) );
            axis square; box off; colormap(b2r(0,1));
            set(gca,'fontname','helvetica','tickdir','out');
            ylabel([num2str(rois.orientations(ori)) '°']);
                        
            subplot_tight(9,4,ori+12,[0.03 0.03]);
            stdshadenew(xvals, mean( Resp{ROI}{ori},1), std( Resp{ROI}{ori},1)./sqrt(6),0.3,'k');
            line([4,4],[0,100], 'linestyle','--','color','r','linewidth',0.5);
            xlim([0,6]);
            ylim([0,max(max(rois.ONResponse{curROIs(ROI)}))]);
            axis square; box off;
            set(gca,'fontname','helvetica','tickdir','out');
            ylabel([num2str(rois.orientations(ori)) '°']);
            
        end;
        
        subplot_tight(9,4,[25 36],[0.05 0.01]);
        errorbar( linspace(0,330,12), mean(rois.shiftON{curROIs(ROI)},2), std(rois.shiftON{curROIs(ROI)},[],2)./sqrt(Params.repTrials),'o','markerfacecolor','k','color','k');
        hold on;
        
        plot( Angles, rois.VMFun{curROIs(ROI)}, 'color','r','linewidth',3); hold on;
        axis square; box off;
        set(gca,'fontname','helvetica','tickdir','out');
        ylabel(['ROI #',num2str(rois.Label(ROI))]);
        set(gca,'xtick', rois.orientations, 'XTickLabel', rois.tcXlabels);
        xlim([-10,360]);
        textYpos = 0.1*ylim;
        textXpos = 0.1*xlim;
        text(9*textXpos(2),5*textYpos(2),['p = ' num2str(rois.sigOriResp(curROIs(ROI)))]);
        text(9*textXpos(2),4*textYpos(2),['GOF = ' num2str(rois.GOF(curROIs(ROI)))]);
        
        clc;
        suptitle([rois.fname(1:end-9),' (ROI# ',num2str(rois.Label(ROI)),')']);
    end;
end

%% Plot Summary TC
if settings.plotSummaryTC
    maxRows = 2;
    subCols = 5;
    maxSubPlots = maxRows * subCols;
    subRows = min([maxRows, ceil(NumROIs/5)]);
    Angles = linspace(0,330, 20000);
    AngRad = Angles.*(pi/180);
    
    pgCt = 1;
    subPlotIndx = 0;
    
    for ROI = 1:NumROIs;
        
        f(pgCt) = figure(200+pgCt); hold on;
        set(gcf,...
            'color','w',...
            'PaperOrientation','landscape',...
            'Position', [10 10 1200 600]);
        subPlotIndx = subPlotIndx +1;
        subplot_tight(subRows,subCols,subPlotIndx,[0.03 0.03]); hold on;
        errorbar( linspace(0,330,12), mean(rois.shiftON{curROIs(ROI)},2), std(rois.shiftON{curROIs(ROI)},[],2)./sqrt(Params.repTrials),'o','markerfacecolor','k','color','k');
        plot( Angles, rois.VMFun{curROIs(ROI)}, 'color','r','linewidth',3); hold on;
        axis square; box off;
        set(gca,'fontname','helvetica','tickdir','out');
        ylabel(['ROI #',num2str(rois.Label(ROI))]);
        set(gca,'xtick', rois.orientations, 'XTickLabel', rois.tcXlabels);
        xlim([-10,360]);
        title(['p = ' num2str(rois.sigOriResp(curROIs(ROI))), ', GOF = ' num2str(rois.GOF(curROIs(ROI)))]);
        
        % Reset subplots
        if rem(ROI,maxSubPlots) == 0
            subPlotIndx = 0;
            pgCt = pgCt + 1;
        end
    end

    % Save figures...
    if settings.saveFigs
        for i = 1:length(f)
            
            suptitle([rois.fname(1:end-9), ' (pg. ', num2str(i), ')']);
            if strcmp(settings.subROIs,'All')
                saveName = [rois.fname(1:end-9) '-allTC-', num2str(i)];
            else
                saveName = [rois.fname(1:end-9) '-filtTC-', num2str(i)];
            end
            saveas(f(i),saveName,'eps2c');
            saveas(f(i),saveName,'jpeg');
            fprintf('Summary TCs figure saved!\n');
            msgCt(3) =  1;
        end
    end
end

%% Plot Filtered ROI TCs
% plotFiltTC = questdlg('Plot filtered ROIs?');

% if strcmp(plotFiltTC,'Yes')
if  settings.filtRoiFlag
    filtQ = inputdlg('Select rois to plot (leave empty for all current ROIs): ', 'Plot individual ROIs', 1, {num2str(rois.filtROI)});
    if isempty(filtQ{:});
        rois.filtROI = rois.Label;
    else
        rois.filtROI = str2num(cell2mat(filtQ));
    end
else
    rois.filtROI = rois.Label;
end

% Plot individual TC for filtered rois
fPlotTC = plotTC(rois, Params, rois.filtROI);

% Save figures...
if settings.saveFigs
    for i = 1:length(fPlotTC)
        saveName = [rois.fname(1:end-9) '-TC-ROI_', num2str(rois.filtROI(i))];
        saveas(fPlotTC(i),saveName,'eps2c');
        saveas(fPlotTC(i),saveName,'jpeg');
    end
    fprintf('Individual TC figure(s) saved! ROIs: %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f \n', rois.filtROI(:))
    fprintf('\n');
    msgCt(4) =  1;
end

%% Plot filtered ROI image.................................................

if settings.plotOverlay
    % Look for jpg of avg projection for movie...
    try
        img = ['AVG_', rois.fname(1:end-9),'.jpg'];
    catch
        fprintf('Select jpg of movie projection...\n');
        img = uigetfile('*.jpg', 'Select projection file (jpg)...');
    end
    
    % Look for zip file containing imagej roi data...
    try
        imageJrois = rois.positions;
    catch
        try
            zipFname = [rois.fname(1:end-8),'RoiSet.zip'];
            imageJrois = ReadImageJROI(zipFname);
        catch
            imageJrois = ReadImageJROI;
        end
    end
    
    %Select appropriate rois to process...
    if settings.filtRoiFlag
        selected = rois.filtROI ;
        selectROIs = imageJrois(selected);
    else
        selected = rois.Label;
        selectROIs = imageJrois(selected);
    end
    
    if settings.plotOSI
        % Plot OSI values for each roi...
        roisVals = rois.OSI(selected);
        maxVal = round(max([0.2,roisVals]),3);
        [f(1)] = labelROIs(img, selectROIs, settings.roiType, roisVals, maxVal);
        colormap('jet');
        myColorbar(maxVal);
        title(rois.fname(1:end-9));
        suptitle(['Orientation Selectivity Index (', settings.subROIs,')']);
        
        % Save figure...
        if settings.saveFigs
            saveName = [rois.fname(1:end-9) '-OSImap (', settings.subROIs,')'];
            saveas(gcf,saveName,'eps2c');
            saveas(gcf,saveName,'jpeg');
            fprintf('OSI rois figure saved!\n');
            msgCt(5) =  1;
        end
    end
    
    if settings.plotPO
        % Plot preferred orientation (PO) values for each roi...
        roisVals = rois.PO(selected);
        maxVal = 180;
        [f(2)] = labelROIs(img, selectROIs, settings.roiType, roisVals, maxVal);
        colormap('jet');
        myColorbar(maxVal,Params.numOris/2);
        title(rois.fname(1:end-9));
        suptitle(['Preferred Orientation (', settings.subROIs,')']);
        
        % Save figure...
        if settings.saveFigs
            saveName = [rois.fname(1:end-9) '-POmap (', settings.subROIs,')'];
            saveas(gcf,saveName,'eps2c');
            saveas(gcf,saveName,'jpeg');
            fprintf('PO rois figure saved!\n');
            msgCt(6) =  1;
        end
    end
end

%% Finish
disp(rois)

if settings.closeAll
    ca;
end
clearvars -except rois