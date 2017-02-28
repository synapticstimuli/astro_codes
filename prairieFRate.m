function [oriFrameRate, newFrameRate] = prairieFRate

%% Attempt to extract frame rate from cfg or env file (PrairieView format.) If  
%  file cannot be read, prompt user for frame rate

directory = dir;

for j = 1:length(directory)
    if ~isempty(strfind(directory(j).name,'cfg')) || ~isempty(strfind(directory(j).name,'env'))
        cfg_idx(j) = 1;
    else
        cfg_idx(j) = 0;
    end
end

% Find the prairie files in the directory
if sum(cfg_idx)==1
    %If only one file exists, set filename to that
    cfg_filename = directory(cfg_idx==1).name;
elseif sum(cfg_idx)>1
    % If multiple files exist, use the first one 
    prairieFileIdxs = find(cfg_idx==1);
    cfg_filename = directory(prairieFileIdxs(1)).name;
    % If multiple files exist, ask user to select correct file
%     cfg_filename = uigetfile({'*.cfg;*.env'},'Open config file');
else
    cfg_filename = [];
end

% Find the frame rate by searching for the correct strinf within the file
% prairie file. Then round to create integer framerate value
if ~isempty(cfg_filename)
    cfg_file = importdata(cfg_filename);
    for j = 1:length(cfg_file)
        if strfind(cfg_file{j},'repetitionPeriod') > 0
            cfg_line = cfg_file{j};
            index = strfind(cfg_line,'repetitionPeriod');
            oriFrameRate = 1/sscanf(cfg_line(index:end),'repetitionPeriod="%f"');
            newFrameRate = ceil(oriFrameRate);
        end
    end
else
    % If no file was found, ask for manual input of new framerate
    oriFrameRate = NaN;
    newFrameRate = input('Could not read file. Enter frame rate (fps): \n');
end
