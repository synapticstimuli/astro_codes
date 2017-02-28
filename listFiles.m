function listFiles(~)


destext = 'First, set parameters for multifile analysis. Then, select folder with experiments. ';
[inputParams, button] = settingsdlg(...
    'Description'                                       , destext                               ,...
    'title'                                             ,'Multifile analysis'                   ,...
    {'Stim Protocol','stim'}                            ,{'Astro','Soma'}                       ,...
    {'Exp type (''NS'',''Spo'')','protocol'}            ,{'NS','Spo'}                           ,...
    {'Channel','channel'}                               ,{'Ch1','Ch2'}                       ,...
    {'Analysis (0=load, 1=rel)','analysis'}             ,{'0','1'}                       ...
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
pause(3)

%% Start 
parfor f = 1:numFolders
    curDir = [folder filesep folderNames{f}];
    cd(curDir)
    curFile = dir([curDir filesep '*' tifFiles]);
    matFiles = dir([curDir filesep '*' inputParams.channel '-reg2_data.mat']);
    
    if inputParams.analysis == 0 
        if isempty(matFiles)
            if ~isempty(curFile)
                fname = curFile.name;
                fprintf('Current file: %s\n',fname);
                fprintf('Creating data file...\n');
                loadMovie_New(inputParams.protocol,fname,curDir);
                pause(0.5);
            else
                fprintf('No registered tif file found in %s, skipping...\n',curDir);
            end
        end
    end
    if inputParams.analysis == 1 
        if ~isempty(matFiles)
%             fprintf('Data file found for %s......skipping.\n',fname);
            relFiles = dir([curDir filesep '*' inputParams.channel '-reg2_RelNat.mat']);
            if isempty(relFiles)
                fprintf('Running rel analysis...\n');
                datafname = strrep(fname,'.tif','_data.mat');
                fprintf('Data file: %s\n', datafname);
%                 diary off
                msg = ['VizReliabilityMap running on:', fname];
                h = msgbox(msg,['Currently running (',num2str(f),' of ',num2str(numFolders),'):']);
                vizReliabilityMap_RG2(datafname);
                close (h)
%                 diary on; disp(datestr(now));
            elseif ~isempty(relFiles)
                fprintf('Reliability already processed for %s...\n',fname);
                if redoFlag
                    vizReliabilityMap_RG2(fname)
                end
            end
        end
    end
end
% diary off
% winopen(diaryFile)
