function imgs_3d_matrix = TIFFreader(tifStackFilePath, format)
    if nargin < 2
        format = 'uint16';
    end
    info = imfinfo(tifStackFilePath);
    numSlices = numel(info);
    firstSlice = imread(tifStackFilePath, 'Index', 1, 'Info', info);
    [height, width] = size(firstSlice);
    imgs_3d_matrix = zeros(height, width, numSlices, format);
    for k = 1:numSlices
        imgs_3d_matrix(:,:,k) = imread(tifStackFilePath, 'Index', k, 'Info', info);
    end
end