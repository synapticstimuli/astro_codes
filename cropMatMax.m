function [croppedMat] = cropMatMax(mat,maxSize)

fprintf('Cropping...\n');
oriSize = size(mat,1);
cropSize = (oriSize - maxSize)/2;
cropSt = cropSize + 1;
cropEnd = cropSize;
matSz = numel(size(mat));

if oriSize == maxSize
    fprintf('No crop needed (size=%1.0f)!\n',oriSize);
    croppedMat = mat;
    return
end

if iscell(mat)
    numCells = length(mat);
    for i = 1:numCells;
        croppedMat{i} = mat{i}(cropSt:end-cropEnd,cropSt:end-cropEnd,:);
    end
elseif ~iscell(mat)
    if matSz==2
        croppedMat = mat(cropSt:end-cropEnd,cropSt:end-cropEnd);
    elseif matSz==3
        croppedMat = mat(cropSt:end-cropEnd,cropSt:end-cropEnd,:);
    end
end

fprintf('Cropped...\n');