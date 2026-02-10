function mask = NucleusMask2(img, close_size)
if nargin<2
    close_size = 3;
end
    if length(size(img)) == 2
        temp = imbinarize(rescale(img), "global");
        temp = imfill(temp, 'holes');
        ccomp = bwconncomp(temp);
        mask = zeros(size(img));
        [~, i] = max(cellfun(@length, ccomp.PixelIdxList));
        mask(ccomp.PixelIdxList{1, i}) = 1;
        mask = imclose(mask, strel('disk',close_size));
    elseif length(size(img)) == 3
        numberOfPages = size(img, 3);
        mask = zeros(size(img));
        for frame_iter = 1:numberOfPages
            temp = img(:, :, frame_iter);
            temp = imbinarize(rescale(temp), "global");
            temp = imfill(temp, 'holes');
            ccomp = bwconncomp(temp);
            temp_mask = zeros(size(temp));
            [~, i] = max(cellfun(@length, ccomp.PixelIdxList));
            if ~isempty(i)
            temp_mask(ccomp.PixelIdxList{1, i}) = 1;
            temp_mask = imclose(temp_mask, strel('disk',close_size));
            end
            mask(:, :, frame_iter) = temp_mask;
        end
    end
    mask = logical(mask);
end