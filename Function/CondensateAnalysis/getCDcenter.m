function CDcenter = getCDcenter(img, img_class, resize_scale, sep_point)
    % use the image before rescaling to get the stable local maxima
    CDcenter = imregionalmax(img, 8);
    CDcenter = imresize(CDcenter, resize_scale);
    CDcenter = CDcenter & img_class>sep_point;
end

