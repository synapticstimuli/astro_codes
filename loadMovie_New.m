function [rois] = loadMovie_New(protocol,fname,pathname)
% Run this script to load and extract the raw data from tiff files.
% Generates the data structure, which includes the dFF, and the movie.dat
% file.
% 
%% Input arguments
rois.WindowSize =   2;          % n X n pixels to bin image
saveFlag        =   1;          % Save data mat file
rois.stimProtocol = 'Astro';
movieFlag = 0;
% if nargin > 3 
%     movieFlag = 1;
% end

if nargin==0
    [fname, pathname] = uigetfile('.tif','Select .tif... ');
    cd(pathname);
    if ismember ('NS',fname)            % Nat Scene v Gratings file
        protocol        =   'NS';
    end
end
cd(pathname);
switch rois.stimProtocol
    case 'Astro'
        rois.stimDur    =   984;
    case 'Soma'
        rois.stimDur    =   792;
end
 
fprintf('Stim protocol: %s-%s...\n', protocol, rois.stimProtocol);
pause(1)

%% Determine frame information from tif file.

try
    info = imfinfo (fname);
    rois.numFrames = numel(info);
    rois.filename = fname;
    rois.yPixels=info(1).Height;
    rois.xPixels=info(1).Width;
catch
    rois.numFrames = input('Enter number of frames: ');
    rois.xPixels = input('Enter number of pixels (X dimension): ');
    rois.yPixels = input('Enter number of pixels (Y dimension): ');
end

%% Attempt to extract frame rate from cfg file (PrairieView format)
    % If cfg file cannot be read, prompt user for frame rate
    directory = dir;
    
    % Determine frame rates...
    [rois.rawFrameRate, rois.newFrameRate] = prairieFRate;

%% Create movie...
if ~movieFlag
    fprintf('Creating movie...\n');
    tic
    tempTif = Tiff(fname, 'r');
    movie = zeros(rois.xPixels,rois.yPixels,rois.numFrames,'uint16');
    
    for i=1:rois.numFrames
        clc
        fprintf('Creating movie...\n');
        fprintf('Running Frame: %d of %d (%0.2fs)%\n',i,rois.numFrames, toc);
        tempTif.setDirectory(i);
        movie(:,:,i)=tempTif.read();
    end
    tempTif.close();
    
    dur(1) = toc/60;
    fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(1));
    pause(1)
elseif movieFlag
    fprintf('Loading movie...\n');
    try
        movieName = strrep(fname, '.tif', '_movie.mat');
    catch
        fprintf('Could not find movie.mat file. Pleae load...\n');
        movieName = uigetfile('*movie.mat');
    end
    load(movieName);
    foo = rem(rois.xPixels,rois.WindowSize);
    if foo~=0
        fprintf('Cropping movie...\n');
        fooEnd = rois.xPixels-foo/2;
        croppedMovie = movie;
    	movie = croppedMovie(foo:fooEnd,foo:fooEnd,:);
        rois.xPixels = size(movie,1);
        rois.yPixels = size(movie,2);
    end
end


 %% (PREVIOUS VER) Interpolate frame rate and calculating dFF
fprintf('Interpolating and obtaining DFF...\n');
tic
oldVector = linspace(0, rois.stimDur, rois.numFrames);
newVector = linspace(0, rois.stimDur, rois.stimDur * rois.newFrameRate);

NumBins = rois.xPixels/rois.WindowSize;
ct = 0;
rois.dFF = zeros(NumBins,NumBins,length(newVector));
for i = 1:NumBins
    for j = 1:NumBins
        ct = ct+1;
        clc
        fprintf('Interpolating and obtaining DFF...\n');
        fprintf('Completed:%0.2f%% (%0.2fs)%\n ',(ct/(NumBins*NumBins))*100,toc);
        foo = movie(1+(i-1)*rois.WindowSize:i*rois.WindowSize,1+(j-1)*rois.WindowSize:j*rois.WindowSize,:);
        foo = mean(mean(foo,1),2);
        foo = spline(oldVector,squeeze(foo),newVector);
        
        final_dFF = computeDFF(foo,rois.newFrameRate);
        rois.dFF(i,j,:) =  final_dFF;
        clear KSD Xi F0 foo
    end;
end;

dur(2) = (toc/60);
fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(2)');
pause(2);

%% Create roi grid positions
if movieFlag
    rois.positions = makeGridMask(rois.xPixels, rois.WindowSize);
end
%% Save data
disp(rois)
switch movieFlag
    case 0
        data = rois;
        fprintf('Saving...\n');
        fprintf('Saving data...\n');
        fname = [pathname filesep rois.filename(1:end-4), '_data'];
        save (fname, 'data', 'rois',  '-v7.3')
        fprintf('Saving movie...\n');
        fname = [pathname filesep rois.filename(1:end-4), '_movie'];
        save (fname, 'movie', '-v7.3')
    case 1
        fprintf('Saving...\n');
        fprintf('Saving data...\n');
        saveName = strrep(rois.filename,'.tif', ['-sw',num2str(rois.WindowSize), '-rois']);
        save ([pathname filesep saveName], 'rois',  '-v7.3')
end

dur(3)=dur(1)+dur(2);
fprintf('Total elapsed time: %0.2fmins\n',dur(3));

