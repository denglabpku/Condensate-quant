function [rna_bkg,base_bkg] = IntensityCalculation(img_series,nucleus_mask, rc_index, rna_width, bkg_width)

numberOfPages = size(img_series, 3);

if nargin < 4
    rna_width = 2; bkg_width = 3;
end
base_bkg = zeros(1, numberOfPages);
rna_bkg = zeros(numberOfPages, 2);

for frame_iter = 1:numberOfPages
    temp_img = double(imgaussfilt(img_series(:, :, frame_iter), 1));
    temp_nucleus_mask = nucleus_mask(:, :, frame_iter);
    % base_bkg(frame_iter) = mean(temp_img(temp_nucleus_mask));
    base_bkg(frame_iter) = prctile(temp_img(:), 90);
    rna_mask = logical(zeros(size(temp_img)));
    rna_mask(ceil(rc_index(frame_iter, 1)), ceil(rc_index(frame_iter, 2))) = 1;
    bkg_mask = imdilate(rna_mask, true(2*bkg_width+1, 2*bkg_width+1));
    rna_mask = imdilate(rna_mask, true(2*rna_width+1, 2*rna_width+1));
    bkg_mask = bkg_mask&~rna_mask;

    bkg = mean(temp_img(bkg_mask));
    rna = sum(temp_img(rna_mask)-bkg);
    rna_bkg(frame_iter, :) = [rna, bkg];
end

end