function croppedMat = cropMat(mat,cropSize)

fprintf('Cropping...\n');
cropSt = cropSize + 1;
cropEnd = cropSize;
matSz = numel(size(mat));

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