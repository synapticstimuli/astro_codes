function roiAnalysis_exp(fname)

fprintf('rois file to load: %s',fname);
load(fname);
fprintf('Loading rois file...\n');

%% Begin ui params input

settings.stimProtocol = 'Astro';
settings.repTrials = 6;
settings.repMovies = 4;
settings.numOris = 12;
settings.runRelAnalysis = 0;
settings.runOriAnalysis = 1;
settings.shiftFlag = 1;
settings.plotTableFlag = 1;
settings.proROI = 1;
settings.subROIs = 'All';
settings.filtRoiFlag = 1;
settings.p = 0.0500;
settings.gofThresh = 0.5000;
settings.osiThresh = 0.1000;
settings.saveROIflag = 1;
settings.saveFigs = 0;
settings.clear = 0;
settings.runROIflag = 0;

% [settings, button] = settingsdlg(...
%     'Description'                                       ,'Set parameters:'                      ,...
%     'title'                                             ,'ROI NM and Ori analysis options'      ,...
%     'separator'                                         ,'Set experiment params:'               ,...
%     {'Select Protocol','stimProtocol'}                  ,{'Astro','Soma'}                       ,...
%     {'Number of trials','repTrials'}                    ,'6'                                     ,...
%     {'Number of movie reps','repMovies'}                ,'4'                                     ,...
%     {'Number of orientations','numOris'}                ,'12'                                    ,...
%     'separator'                                         ,'Set analysis params:'                 ,...
%     {'Calculate NM reliability?','runRelAnalysis'}      ,[true]                                 ,...
%     {'Compute Ori tuning curves?','runOriAnalysis'}     ,[true]                                 ,...
%     {'Shift TC axes to center at 0°?','shiftFlag' }     ,[true]                                 ,...
%     {'Plot table?', 'plotTableFlag'}                    ,[true]                                 ,...
%     'separator'                                         ,'Select ROIs:'                         ,...
%     {'Process converted rois?', 'proROI'}               ,[false]                          ,...
%     {'Select ROIs to run','subROIs'}                   ,'All'                                   ,...
%     'separator'                                         ,'Set filtering params:'               ,...
%     {'Filter ROIs?','filtRoiFlag'}                      ,[false]                                ,...
%     {'Response p value threshold','p'}                 ,'0.05'                                  ,...
%     {'GOF threshold','gofThresh'}                      ,'0.5'                                   ,...
%     {'OSI threshold','osiThresh'}                      ,'0.1'                                   ,...
%     'separator'                                         ,'Save and Finish:'                     ,...
%     {'Save rois.mat?','saveROIflag'}                    ,[true]                                 ,...
%     {'Save all figures?','saveFigs'}                    ,[false]                                 ,...
%     {'Clear all?', 'clear'}                             ,[false]                                ,...
%     {'Run roisPlot?','runROIflag'}                      ,[false]                                 ...
%     );

% if strcmp(button,'cancel') || isempty(button)
%     fprintf('\n*****************  No settings were entered. Script stopped.  *********************\n\n\n');
%     return
% end
Params.stimProtocol     = settings.stimProtocol;
Params.repTrials        = settings.repTrials;
Params.repMovies        = settings.repMovies;
Params.numOris          = settings.numOris;
%% Load Data...............................................................
if ~settings.proROI
    % Load roi raw data...
%     rois = [];
    [rois, NumROIs, rois.fname] = loadROI;
    rois.Label = 1:NumROIs;
    
    % Determine frame rates...
    [rois.rawFrameRate, rois.newFrameRate] = prairieFRate;
else
    % Select and load converted rois file.........
    if ~exist('rois','var')
        roisFile = uigetfile('*rois.mat','Select rois.mat file...');
        
        if roisFile==0
            fprintf('No rois.mat file selected. Terminating script.\n');
            return
        else
            fprintf('Loading rois.mat...\n');
            load(roisFile,'rois');
            disp(rois);
        end
        rois.fname = roisFile;
    end
    rois.filename = fname;
    NumROIs = size(rois.dFFNew,1);
    rois.Label = 1:NumROIs;
end

% Modify roi labels and indices for pre-selected rois
if  ~strcmp(settings.subROIs,'All')
    subROIs     = str2num(settings.subROIs);
    rois.rawF   = rois.rawF(subROIs,:);
    NumROIs     = size(rois.rawF,1);
    rois.Label  = rois.Label(subROIs);
    settings.filtRoiFlag = 0;
    settings.saveROIflag = 0;
end

%% Stimulus parameters and Housekeeping
switch Params.stimProtocol
    case 'Astro'
        rois.stimDur           = 984;
        Params.durBlank_nm     = 4;
        Params.durBlank_grats  = 6;
        Params.durStim_ns      = 13;
        Params.durStim_grats   = 2;
    case 'Soma'
        rois.stimDur           = 792;
        Params.durBlank_nm     = 2;
        Params.durBlank_grats  = 4;
        Params.durStim_ns      = 13;
        Params.durStim_grats   = 2;
end
rois.orientations = 0 : 360/Params.numOris :330;     % Range of degrees based on numOris

%% Compute dFF.............................................................

if ~settings.proROI
    fprintf('Computing DFF...\n');
    for i = 1:NumROIs
        rois.dFF(i,:) = computeDFF(rois.rawF(i,:),rois.newFrameRate);
    end;
    fprintf('...Done!\n');
end

%% Interplolate to new frame rate
if ~settings.proROI
    tic;
    ct = 0;
    
    NumSamplsNew = rois.stimDur .* rois.newFrameRate;
    OldVec       = linspace(0, rois.stimDur, size(rois.dFF,2));
    NewVec       = linspace(0, rois.stimDur, NumSamplsNew);
    
    for i = 1:NumROIs
        % Display progress...
        ct = ct+1;
        clc
        fprintf('Interpolating DFF...\n');
        fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/NumROIs)*100,toc);
        
        % Create new dff vectors...
        dFFNew(i,:) =  spline(OldVec, rois.dFF(i,:), NewVec);
    end;
    
    rois.dFFNew = dFFNew;
    clearvars OldVec NewVec ct
    fprintf('...Done!\n');
end

%% Break into Trials
tic
ct = 0;

% Trial specific structures...
NumSamplsNew = rois.stimDur .* rois.newFrameRate;                                       % Total number of frames for entire experiment
TimeEachTrial   = (Params.repMovies * (Params.durStim_ns + Params.durBlank_nm))...      % Duration (s) of each trial
                  + (Params.numOris * (Params.durStim_grats + Params.durBlank_grats));
SamplsEachTrial = TimeEachTrial * rois.newFrameRate;                                    % Number of frames for each trial
NumTrials       = NumSamplsNew ./ SamplsEachTrial;                                      % Calculation of num trials (see also Params.repTrials)
Trials          = reshape( rois.dFFNew, [size(rois.dFFNew,1), SamplsEachTrial, NumTrials] );      % Reshape dff into frames of trial * num trials

% Stimulus-specific trial structures... 
durSec_NM       = Params.durStim_ns + Params.durBlank_nm;                                   % Duration (s) of NS (ON + OFF)
durFrames_NM    = durSec_NM * rois.newFrameRate;                                       % Num frames of NS stim (ON + OFF)
durTrials_NM    = durFrames_NM * Params.repMovies;                                     % Total length of NS (ON + OFF) per trial
durSec_Grats    = Params.durStim_grats + Params.durBlank_grats;                        % Duration (s) of Grats (ON + OFF)
durFrames_Grats = durSec_Grats * rois.newFrameRate;                                 % Num frames of Grats stim (ON + OFF)
durTrials_Grats = durFrames_Grats * Params.numOris;
numTrials_Grats = Params.numOris * Params.repTrials;                                % Total numer of gratings repetitions

NatMovies = [];
Gratings = [];

for i = 1:NumTrials
    % Display progress...
    ct = ct+1;
    clc
    fprintf('Partitioning data...\n');
    fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/NumTrials)*100,toc);
    
    % Partitioning...
    foo = squeeze( Trials(:,:,i) );
    NatMovies(:,:,i) =  foo(:, 1:durTrials_NM);
    Gratings(:,:,i)  =  foo(:, durTrials_NM+1:end);
end;

% Store values
rois.NatMovies = NatMovies;
rois.Gratings = Gratings;

clearvars ct button fo Gratings NatMovies roisFile

%% Calculate NM reliabiity................................................
if settings.runRelAnalysis
    if ~settings.proROI
        tic;
        ct = 0;
        
        respTest = [];
        Reliability = NaN(1,NumROIs);
        
        %Set permutation
        numTrials_NM = Params.repMovies * NumTrials;
        P  = nchoosek( 1:numTrials_NM, 2);
        
        for ROI = 1:NumROIs;
            % Display progress...
            ct = ct+1;
            clc
            fprintf('Computing NM Reliability....\n');
            fprintf('Completed:%0.3f%% (%0.2f min)\n',(ct/NumROIs)*100,toc/60);
            
            resp_NM = squeeze( rois.NatMovies(ROI,:,:) );
            resp_NM = resp_NM(:);
            respTr_NM = reshape( resp_NM', [durFrames_NM, Params.repTrials.*Params.repMovies] );
            respTr_NM_norm = respTr_NM./max(respTr_NM(:));
            
            template = mean(respTr_NM,2)';
            respTest(ROI,:) = template(1,:);
            
            % Calculate reliability accross trials
            for j = 1:size(P,1);
                T1 = P(j,1);
                T2 = P(j,2);
                R(j) =  distcorr( respTr_NM(:,T1), respTr_NM(:,T2) );
            end;
            
            Reliability(ROI) = mean(R);
            
            %     for j = 1:24
            %         Rel(j) = corr( respTr_NM(:,j), template');
            %     end;
            
        end;
        
        % Store values
        rois.Reliability = Reliability;
        
        clearvars P ct
    end
    
end

%% Compute Ori TC from Grating.....
if settings.runOriAnalysis
    progressbar;
    for ROI = 1:NumROIs;
        resp_Grat = squeeze( rois.Gratings(ROI,:,:) );
        resp_Grat = resp_Grat(:);
        respTr_Grat     = reshape( resp_Grat', [durFrames_Grats, numTrials_Grats] );
        respTr_Grat_norm = respTr_Grat./max(respTr_Grat(:));
        
        for ori = 1:Params.numOris
            ct = 0;
            for ind = ori : Params.numOris : 72
                ct = ct+1;
                Resp{ROI}{ori}(ct,:) = respTr_Grat( :, ind)';
            end;
            
            ON{ROI}(ori,:) = ( mean( Resp{ROI}{ori}(:, ((durSec_Grats-Params.durStim_grats)*rois.newFrameRate):(durFrames_Grats)), 2) )';
            Tun{ROI}(ori,:) = ON{ROI}(ori,:);
            OFF{ROI}(ori,:)  = mean(  Resp{ROI}{ori}(:, (2) * rois.newFrameRate:(durSec_Grats-Params.durStim_grats)*rois.newFrameRate), 2)';
            STDOFF{ROI}(ori) = mean(std( Resp{ROI}{ori}(:, (2) * rois.newFrameRate:(durSec_Grats-Params.durStim_grats)*rois.newFrameRate),0,2));
        end;
        
        if settings.shiftFlag
            zeroCent = [8 9 10 11 12 1 2 3 4 5 6 7]';
            currON = [zeroCent ON{ROI}];
            currSorted = sortrows(currON);
            shiftON{ROI} = currSorted(:,2:end);
            tcXlabels = {'','','-90','','','0','','','90','','','180'};
            
        else
            shiftON{ROI}=ON{ROI};
            tcXlabels = {'0','','','90','','','180','','','270','',''};
        end
        rois.shiftON = shiftON;
        rois.tcXlabels = tcXlabels;
        
        % Soma.TCFit   = F;
        rois.TCraw{ROI}   = (mean(ON{ROI},2)'- mean(OFF{ROI},2)');
        rois.TCraw{ROI}   = mean(ON{ROI},2)';
        rois.ONResponse  = ON;
        rois.OFFResponse = OFF;
        [rois.OSI(ROI), rois.PO(ROI), rois.DSI(ROI)] = CalculateOSI(rois.TCraw{ROI});
        rois.maxOriResp(ROI) = max((mean(rois.ONResponse{ROI},2)));
        rois.meanOriResp(ROI) = mean((mean(rois.ONResponse{ROI},2)));
        rois.normOriResp(ROI) = (max(mean(rois.ONResponse{ROI},2))-(mean(mean(rois.OFFResponse{ROI},2))));
        [h p] = ttest(rois.ONResponse{ROI}(:), rois.OFFResponse{ROI}(:));
        rois.sigOriResp(ROI) = p;
        %     for i=1:NumROIs
        %         Filt{ROI}=mean(Tun{ROI}(:));
        %     end
        
        % Run VM function for TC fit
        
        Angles = linspace(0,330, 20000);
        AngRad = Angles.*(pi/180);
        [CoeffSet, GOF] = VMFit(mean(shiftON{ROI},2)',30);
        rois.VMFun{ROI} = VonMisesFuntion(CoeffSet,AngRad);
        progressbar(ROI/NumROIs);
        h_pb = gcf;
        rois.GOF(ROI) = GOF;
        
        clc;
    end;
    close(h_pb); 
    
    % Convert all PO to 180°
    poLarge = find(rois.PO>180);
    rois.PO(poLarge) = rois.PO(poLarge)-180;
end

%% Create table with slope coefficients
if settings.plotTableFlag
    table = figure; hold on;
    saveTable = 1;
    tableData = [];
    
    % Build table...
    for i = 1:NumROIs;
        tableData(i,1) = rois.Label(i);
        tableData(i,2) = rois.Reliability(i);
        tableData(i,3) = rois.OSI(i);
        tableData(i,4) = rois.PO(i);
        tableData(i,5) = rois.DSI(i);
        tableData(i,6) = rois.sigOriResp(i);
        tableData(i,7) = rois.GOF(i);
    end
    supTitleText = [rois.filename(1:end-9),' (All ROIs)'];
    if settings.filtRoiFlag
        % Select filter parameters and filter rois....
        filtROI = find(rois.sigOriResp<settings.p & rois.GOF>settings.gofThresh & rois.OSI>settings.osiThresh);
        numFilt = length(filtROI);
        temp = NaN(numFilt+1,7);
        tableData = [tableData; temp];
        k=0;
        % Build portion of table with just filetered rois...
        for j = NumROIs+2:size(tableData,1);
            k=k+1;
            tableData(j,1) = rois.Label(filtROI(k));
            tableData(j,2) = rois.Reliability(filtROI(k));
            tableData(j,3) = rois.OSI(filtROI(k));
            tableData(j,4) = rois.PO(filtROI(k));
            tableData(j,5) = rois.DSI(filtROI(k));
            tableData(j,6) = rois.sigOriResp(filtROI(k));
            tableData(j,7) = rois.GOF(filtROI(k));
        end
        supTitleText = [rois.filename(1:end-9), ' (p<', num2str(settings.p),...
            ', GOF>', num2str(settings.gofThresh), ', OSI>', num2str(settings.osiThresh), ')'];
        rois.filtROI = filtROI;
    end
    
    % Set column titles and formats...
    colNames = {'ROI' 'Reliability' 'OSI' 'PO' 'DSI' 'Sig (p)' 'GOF'};
    colForm = { 'numeric', 'numeric', 'numeric','numeric','numeric','numeric','numeric'};
    
    % Generate table...
    t = uitable(table,'Data', tableData, 'ColumnName', colNames, 'ColumnFormat', colForm);
    
    % Set width and height...
    t.Position(3) = 560;
    t.Position(4) = 500;
    t.ColumnWidth={'auto','auto','auto','auto','auto','auto','auto'};
    %Set table figure properties...     
    set(table,...
        'color'         ,'w'                                        ,...
        'Position'      ,[100 100 t.Position(3)+50 t.Position(4)+70],...
        'NumberTitle'   ,'off'                                      );

    suptitle(supTitleText);
    
    % Save table as matlab figure...
    tableName = [rois.filename(1:end-9) '-fTable'];
    saveas(gcf,tableName,'fig');
%     saveas(gcf,tableName,'eps2c');
    fprintf('Table (.fig) saved!\n');
end

 %% Save

if settings.saveROIflag
    fprintf('Saving...\n');
    fprintf('Saving...rois...\n');
    saveName = [rois.filename(1:end-4),'.mat'];
    save(saveName, 'Params','rois', '-v7.3');
    fprintf('Saving...rois...Done!\n');
end

%% Run plot

if settings.runROIflag
    
   fprintf('Running plotting script...\n');
   pause(2);
   roiPlot(saveName);
   clearvars -except rois Params
end

%% Clear

if settings.clear
    disp(rois)
    clearvars 
end
