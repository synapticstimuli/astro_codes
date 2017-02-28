function [roisStruct] = makeGridMask(imgSize, gridSize)

if nargin == 0
    fprintf('Select jpg of movie projection...\n');
    imgFile = uigetfile('*.jpg', 'Select projection file (jpg)...');
    img = imread(imgFile);
    imgSize = size(img);
    gridSize = inputdlg('Enter window sze:','Window bin size',1,{'2'});
elseif nargin ==1
    gridSize = 2;
end
if iscell(gridSize)
    gridSize = str2double(gridSize{:});
end
    
sizeX = imgSize(1);
sizeY = imgSize(2);

rows = floor(sizeY/gridSize);
cols = floor(sizeX/gridSize);
numGrid = rows * cols;

xpos = 1:gridSize:sizeX;
ypos = 1:gridSize:sizeY;

roisStruct = cell(1,numGrid);
for s = 1:length(roisStruct)
    roisStruct{s} = struct('vnRectBounds',[],'cols',[],'rows',[]);
    roisStruct{s}.vnRectBounds = zeros(numGrid, 4);
    roisStruct{s}.cols = zeros(numGrid, 4);
    roisStruct{s}.rows = zeros(numGrid, 4);
end

ct = 0;
tic
for i = 1:rows
    for j = 1:cols
        ct = ct +1;
        clc
        fprintf('Creating positons matrix...\n');
        fprintf('Completed: %0.2f%% (%0.2fs)%\n ',(ct/(rows*cols))*100,toc);
        curVectorBounds = [xpos(i) ypos(j) (xpos(i) + gridSize) (ypos(j) + gridSize)];
        curCols = [xpos(i) (xpos(i) + gridSize) (xpos(i) + gridSize) xpos(i)];
        curRows = [ypos(j) ypos(j) (ypos(j) + gridSize) (ypos(j) + gridSize)];
        roisStruct{ct}.vnRectBounds = curVectorBounds;
        roisStruct{ct}.cols = curCols;
        roisStruct{ct}.rows = curRows;
    end
end

