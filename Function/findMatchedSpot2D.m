function matched_spot = findMatchedSpot2D(img, nucleus_mask, spots_405_3D)

    [spots, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    if isempty(spots)
        error('NucleusMask:No405Spot', 'No 405 spot was detected inside the nucleus mask.');
    end
    
    id = abs(spots(:, 1)-spots_405_3D(1))<3 & abs(spots(:, 2)-spots_405_3D(2))<3;
    spots = spots(id, :);
    quality = quality(id);
    if isempty(spots)
        error('NucleusMask:NoMatched405Spot', 'No 2D 405 spot matched the 3D 405 spot inside the nucleus mask.');
    end
    
    [~, i] = max(quality);
    matched_spot = spots(i, :);
    
end