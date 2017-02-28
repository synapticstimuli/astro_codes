 function multiFolderAnalysis(~)

destext = 'First, set parameters for multifile analysis. Then, select folder with experiments. ';
[inputParams, button] = settingsdlg(...
    'Description'                                       , destext                               ,...
    'title'                                             ,'Multifile analysis'                   ,...
    {'Stim Protocol','stim'}                            ,{'Astro','Soma'}                       ,...
    {'Exp type (''NS'',''Spo'')','protocol'}            ,{'NS-','Spo-'}                           ,...
    {'Channel','channel'}                               ,{'Ch2','Ch1'}                          ,...
    {'Analysis','analysis'}                             ,{'loadMovie','vizRelMap','cCorr'       ,...
    'convertROI','roiAnalysis','roiPlot'}                                                                               ...
    );
if strcmp(button,'cancel') || isempty(button)
    fprintf('\n***  No settings were entered. Script stopped.  ***\n\n\n');
    return
end

tifFiles = [inputParams.channel '-reg2.tif']; 
folder = uigetdir;
if ~folder
    fprintf('*** No folder selected. Script aborted ***\n');
end
cd(folder);

%% 

% Find files in the current directory matching the protocol id, eg '-NS-'
clc
folders = dir([folder filesep '*' inputParams.protocol '*' ]);
numFolders = length(folders);

fprintf('%2.0f folders found:\n',numFolders);

% Print found file names
for f = 1: numFolders
    folderNames{f,1} = folders(f).name;
    fprintf('%s\n',folderNames{f,1});
end
pause(2)

%% Start slected analysis

switch inputParams.analysis
    case 'loadMovie'
        parfor f = 1:numFolders
            % Select correct files in each folder...
            curDir = [folder filesep folderNames{f}];
            cd(curDir)
            curFile     = dir([curDir filesep '*' tifFiles]);
            matFiles    = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
            
            if isempty(matFiles)
                if ~isempty(curFile)
                    fname       = curFile.name;
                    fprintf('Current file: %s\n',fname);
                    fprintf('Creating data file...\n');
                    loadMovie_New(inputParams.protocol,fname,curDir);
                else
                    fprintf('No registered tif file found in %s, skipping...\n',curDir);
                end
            end
        end

    case 'vizRelMap'
        parfor f = 1:numFolders
            curDir = [folder filesep folderNames{f}];
            cd(curDir)
            curFile = dir([curDir filesep '*' tifFiles]);
            matFiles = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
            fname = curFile.name;
            if inputParams.analysis == 1
                if ~isempty(matFiles)
                    relFiles = dir([curDir filesep '*' inputParams.channel '-reg2_RelNat.mat']);
                    if isempty(relFiles)
                        fprintf('Running rel analysis...\n');
                        datafname = strrep(fname,'.tif','_data.mat');
                        fprintf('Data file: %s\n', datafname);
                        msg = ['VizReliabilityMap running on:', fname];
                        h = msgbox(msg,['Currently running (',num2str(f),' of ',num2str(numFolders),'):']);
                        vizReliabilityMap_RG2(datafname);
                        close (h)
                    elseif ~isempty(relFiles)
                        fprintf('Reliability already processed for %s...\n',fname);
                        if redoFlag
                            vizReliabilityMap_RG2(fname)
                        end
                    end
                end
            end
            
        end
        
    case 'cCorr'
        parfor f = 1:numFolders
            curDir = [folder filesep folderNames{f}];
            cd(curDir)
            curFile = dir([curDir filesep '*' tifFiles]);
            matFiles = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
            fname = curFile.name;
            if inputParams.analysis == 1
                if ~isempty(matFiles)
                    fprintf('Running cross correlation...\n');
                    datafname = strrep(fname,'.tif','_data.mat');
                    fprintf('Data file: %s\n', datafname);
                    msg = ['CrossCorrImage running on:', fname];
                    h = msgbox(msg,['Currently running (',num2str(f),' of ',num2str(numFolders),'):']);
                    [ccImg,thImg]= CrossCorrImage(threshVal,movie)
                    close (h)
                end
            end
        end
        
    case 'convertROI'
%         progressbar;
%         for f = 1:numFolders
        parfor f = 1:numFolders
            curDir = [folder filesep folderNames{f}];
            cd(curDir)
            curFile = dir([curDir filesep '*' tifFiles]);
            matFiles = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
            fname = curFile.name;
            if ~isempty(matFiles)
                fprintf('Running convert2ROIs...\n');
                datafname = strrep(fname,'.tif','_data.mat');
                fprintf('Data file: %s\n', datafname);
%                 msg = ['convert2ROIs running on:', fname];
%                 h = msgbox(msg,['Currently running (',num2str(f),' of ',num2str(numFolders),'):']);
%                 progressbar(f/numFolders);
%                 h_pb = gcf;
                convert2ROIs(datafname)
            elseif isempty(matFiles)
                fprintf('No reg2_data.mat file found in %s...\n', curDir);
            end
        end
%         close(h_pb);

    case 'roiAnalysis'
        %         progressbar;
        %         for f = 1:numFolders
%         parloops = numFolders;
%         hbar = parfor_progressbar(parloops,'Processing roiAnalysis_exp...');
        parfor f = 1:numFolders
            curDir = [folder filesep folderNames{f}];
            cd(curDir)
            curFile = dir([curDir filesep '*' tifFiles]);
            matFiles = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
            fname = curFile.name;
            if ~isempty(matFiles)
                fprintf('Running roiAnalysis_exp...\n');
                roisfname = strrep(fname,'.tif','-rois.mat');
                fprintf('Data file: %s\n', roisfname);
                %                 msg = ['convert2ROIs running on:', fname];
                %                 h = msgbox(msg,['Currently running (',num2str(f),' of ',num2str(numFolders),'):']);
                %                 progressbar(f/numFolders);
                %                 h_pb = gcf;
                roiAnalysis_exp(roisfname)
            elseif isempty(matFiles)
                fprintf('No rois.mat file found in %s...\n', curDir);
            end
            
        end
        
        %         close(h_pb);
end

