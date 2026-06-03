function CDinterface = getCDinterface(img, img_class, sep_point)
    % interface
    CDinterface = logical(zeros(size(img)));
    temp_img_class = zeros(size(img_class, 1)+2, size(img_class, 2)+2);
    temp_img_class(2:(end-1), 2:(end-1)) = img_class;
    for x = 2:(size(temp_img_class, 2)-1)
        for y = 2:(size(temp_img_class, 1)-1)
            if temp_img_class(y, x) == sep_point+1
                temp_neighbor = temp_img_class((y-1):(y+1), (x-1):(x+1));
                if sum(temp_neighbor(:) == sep_point) > 0
                    CDinterface(y-1, x-1) = 1;
                end
            end
        end
    end
end

