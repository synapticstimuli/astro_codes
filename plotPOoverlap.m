function plotPOoverlap(~)
%% Set parameters.........................................................
stdVal = 5;
imgAlpha = 0.3;
saveFlag = 1;

%% Select rois.mat files..................................................
files = dir('*rois.mat');
if numel(files)> 2
    fprintf('Too many reg2.tif files found.\n');
    files = uigetfile('*rois.mat','Select rois.mat files for each channel...',...
        'MultiSelect','on');
    ch1fname = files{1};
    ch2fname = files{2};
else
    ch1fname = files(1).name;
    ch2fname = files(2).name;
end

%% Load roi files.........................................................
fprintf('Loading files...');
ch1 = load(ch1fname,'rois');
ch2 = load(ch2fname,'rois');
fprintf('done!\n');

 %% Look for jpg of avg projection for movie...
    try
        ch1imgName = ['AVG_', ch1.rois.fname(1:end-9),'.jpg'];
        ch2imgName = ['AVG_', ch2.rois.fname(1:end-9),'.jpg'];
        ch1img = imread(ch1imgName);
        ch2img = imread(ch2imgName);        
    catch
        fprintf('Select jpg of movie projection for ch1...\n');
        ch1imgName = uigetfile('*.jpg', 'Select ch1 projection file (jpg)...');
        fprintf('Select jpg of movie projection for ch2...\n');
        ch2imgName = uigetfile('*.jpg', 'Select ch2 projection file (jpg)...');
    end
    
    ch1img = imresize(ch1img,0.25,'nearest');
    ch2img = imresize(ch2img,0.25,'nearest');
    ch1rgb = cat(3,ch1img,ch1img,ch1img);
    ch2rgb = cat(3,ch2img,ch2img,ch2img);
    
%% Load values for each channel...........................................

ch1vals = ch1.rois.PO;
ch2vals = ch2.rois.PO;
[ch1POidx] = sigPOresp(ch1, stdVal);
[ch2POidx] = sigPOresp(ch2, stdVal);
ch1POidx(~ch1POidx) = NaN;
ch2POidx(~ch2POidx) = NaN;

%% Create colocalization matrix
matSize = ch1.rois.xPixels / ch1.rois.windowSize;
expName = ch1fname(1:strfind(ch1fname,'-Ch1')-1);

% Discretize angles to conform to the number and degree used in th vis
% stim.
values = [];
edges = ch1.rois.orientations;
for i = 2: length(edges)
    values{i-1} = num2str(edges(i));
end
ch1vals_disc = str2num(char(discretize(ch1vals,edges,'categorical',values)'));
ch2vals_disc = str2num(char(discretize(ch2vals,edges,'categorical',values)'));
ch1Idx = ch1POidx;
ch2Idx = ch2POidx;
% ch1Idx = (ch1vals_disc.*ch1POidx') / 180;
% ch2Idx = (ch2vals_disc.*ch2POidx') / 180;



%% Create matrix to plot

ch1Mat = reshape(ch1Idx,matSize,[])';
ch1Matrgb = zeros(matSize,matSize,3);
ch1Matrgb(:,:,1) = ch1Mat;
ch2Mat = reshape(ch2Idx,matSize,[])';
ch2Matrgb = zeros(matSize,matSize,3);
ch2Matrgb(:,:,2) = ch2Mat;
% Create rgb matrix with thresholded channel data
rgbImage = zeros(matSize, matSize, 3);
rgbImage(:,:,1) = ch1Mat;
rgbImage(:,:,2) = ch2Mat;

%% Plot figure...........................................................
f = figure('color','w'      ,...
           'position', [300 300 900 400]);

h(1) = subplot_tight(1,3,1);
imshow(rgbImage);
title('Ch1(red) + Ch2(grren)');
axis square; box off;

h(2) = subplot_tight(1,3,2);
ch1final = (imgAlpha*uint8(255.*ch1Matrgb)) + ch1img;
imagesc(ch1final);
title('Ch1');
axis square; box off;

h(3) = subplot_tight(1,3,3);
ch2final = (imgAlpha*uint8(255.*ch2Matrgb)) + ch2img;
imagesc(ch2final);
title('Ch2');
%  myColorbar(180,6);
axis square; box off;

suptitle([expName, '-(PO-', num2str(stdVal), 'std)']);
contourset(h(1:3),                            ...
    'XTickLabel','',...
    'YTickLabel','',...
    'XMinorGrid','on',...
    'YMinorGrid','on',...
    'GridColor',[1 1 1],...
    'GridAlpha',0.5,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'MinorGridColor',[1 1 0],...
    'MinorGridLineStyle','-',...
    'MinorGridAlpha',0.15,...
    'Layer','top'...
    );

%% Save figure

if saveFlag
    fprintf('Saving figure...\n');
    saveName = [expName, '-(PO-', num2str(stdVal), 'std)'];
    saveas(gcf, saveName, 'jpeg');
    saveas(gcf, saveName, 'eps2c');
end
%% Nested functions.......................................................
    
    % Identify rois with sig PO responses.....
    function [sigPOidx] = sigPOresp(chMat, stdVal)
        
        numROIs = numel(chMat.rois.PO);
        for ni = 1: numROIs
           onOFF = chMat.rois.ONResponse{ni}-chMat.rois.OFFResponse{ni};
           maxOri = max(mean(onOFF,2));
           minOri = min(mean(onOFF,2));
           meanOri = mean(mean(onOFF),2);
           stdOri = std(mean(onOFF));
           sigPOidx(ni) = maxOri-meanOri > stdVal * stdOri;
        end
        
        sigPOidx = double(sigPOidx);
    end

end