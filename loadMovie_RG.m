function [data] = loadMovie_RG(protocol,fname,pathname)
% Run this script to load and extract the raw data from tiff files.
% Generates the data structure, which includes the dFF, and the movie.dat
% file.
% 
%% Input arguments
data.WindowSize =   2;          % n X n pixels to bin image
saveFlag        =   1;          % Save data mat file
displayFlag     =   0;          % Display movie
data.stimProtocol = 'Astro';

if nargin==0
    [fname, pathname] = uigetfile('.tif','Select .tif... ');
    cd(pathname);
    if ismember ('NS',fname)            % Nat Scene v Gratings file
        protocol        =   'NS';
        switch data.stimProtocol
            case 'Astro'
                data.stimDur    =   984;
            case 'Soma'
                data.stimDur    =   792;
        end
%     elseif ismember ('RF',fname)        % Short receptive field 
%         data.stimDur    =   500;
%         protocol        =   'RF';
%     elseif ismember ('SN',fname)        % Sparse noise file
%         data.stimDur    =   1200;
%         protocol        =   'NS';
%     elseif ismember ('Spo',fname)                           % Spontaneouos activity file
%         data.stimDur    =   189;
%         protocol        =   'Spo';
    end
% elseif nargin==2
%     cd(pathname);
%     switch protocol
%         case 'NS'
%             data.stimDur    =   984;        
%         case 'RF'
%             data.stimDur    =   500;
%         case 'SN'
%             data.stimDur    =   1200;
%         case 'Spo'
%             data.stimDur    =   189;
%     end
elseif nargin == 3
    cd(pathname);
    switch protocol
        case 'NS'
            switch data.stimProtocol
                case 'Astro'
                    data.stimDur    =   984;
                case 'Soma'
                    data.stimDur    =   792;
            end
        case 'RF'
            data.stimDur    =   500;
        case 'SN'
            data.stimDur    =   1200;
        case 'Spo'
            data.stimDur    =   189;
    end
end

fprintf('Stim protocol: %s-%s...\n', protocol, data.stimProtocol);
pause(1)

%% Determine frame information from tif file.
try
    info = imfinfo (fname);
    data.numFrames = numel(info);
    data.filename = fname;
    data.yPixels=info(1).Height;
    data.xPixels=info(1).Width;
catch
    data.numFrames = input('Enter number of frames: ');
    data.xPixels = input('Enter number of pixels (X dimension): ');
    data.yPixels = input('Enter number of pixels (Y dimension): ');
end

%% Attempt to extract frame rate from cfg file (PrairieView format)
    % If cfg file cannot be read, prompt user for frame rate
    directory = dir;
    
    for j = 1:length(directory)
        if ~isempty(strfind(directory(j).name,'cfg')) | ~isempty(strfind(directory(j).name,'env'))
            cfg_idx(j) = 1;
        else
            cfg_idx(j) = 0;
        end
    end
    if sum(cfg_idx)==1
        cfg_filename = directory(cfg_idx==1).name;
    elseif sum(cfg_idx)>1
        cfg_filename = uigetfile({'*.cfg;*.env'},'Open config file');
    else
        cfg_filename = [];
    end
    if ~isempty(cfg_filename)
        cfg_file = importdata(cfg_filename);
        for j = 1:length(cfg_file)
            if strfind(cfg_file{j},'repetitionPeriod') > 0
                cfg_line = cfg_file{j};
                index = strfind(cfg_line,'repetitionPeriod');
                data.rawFrameRate = 1/sscanf(cfg_line(index:end),'repetitionPeriod="%f"');
                data.newFrameRate = ceil(data.rawFrameRate);
            end
        end
    else
        data.newFrameRate = input('Enter frame rate (Hz): ');
    end
    
% displayData(data);
% fprintf('Data struct attributes:\n');
% disp(data)
% pause(3)

%% Create movie...
tic
 
tempTif = Tiff(fname, 'r');
movie = zeros(data.xPixels,data.yPixels,data.numFrames,'uint16');

for i=1:data.numFrames
    clc
    fprintf('Creating movie...\n');
    fprintf('Running Frame: %d of %d (%0.2fs)%\n',i,data.numFrames, toc);
    tempTif.setDirectory(i);
    movie(:,:,i)=tempTif.read();
end
tempTif.close();

dur(1) = toc/60;
fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(1));
pause(2)

%% Interpolate frame rate and calculating dFF
fprintf('Interpolating and obtaining DFF...\n');
tic
oldVector = linspace(0, data.stimDur, data.numFrames);
newVector = linspace(0, data.stimDur, data.stimDur*data.newFrameRate);

NumBins = data.xPixels/data.WindowSize;
ct = 0;
data.dFF = zeros(NumBins,NumBins,length(newVector));
for i = 1:NumBins
    for j = 1:NumBins
        ct = ct+1;
        clc
        fprintf('Interpolating and obtaining DFF...\n');
        fprintf('Completed:%0.2f%% (%0.2fs)%\n ',(ct/(NumBins*NumBins))*100,toc);
        foo = movie(1+(i-1)*data.WindowSize:i*data.WindowSize,1+(j-1)*data.WindowSize:j*data.WindowSize,:);
        foo = mean(mean(foo,1),2);
        foo = spline(oldVector,squeeze(foo),newVector);
        
        final_dFF = computeDFF(foo,data.newFrameRate);
        data.dFF(i,j,:) =  final_dFF;
        clear KSD Xi F0 foo
    end;
end;

dur(2) = (toc/60);
fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(2)');
pause(2);

%% Save data
disp(data)

if saveFlag
    fprintf('Saving...\n');
    fprintf('Saving data...\n');
    fname = [pathname filesep data.filename(1:end-4), '_data'];
    save (fname, 'data',  '-v7.3')
    fprintf('Saving data...Done.\n');
    fprintf('Saving movie...\n');
    fname = [pathname filesep data.filename(1:end-4), '_movie'];
    save (fname, 'movie', '-v7.3')
    fprintf('Saving movie...Done.\n');
else
end
dur(3)=dur(1)+dur(2);
fprintf('Total elapsed time: %0.2fmins\n',dur(3));

