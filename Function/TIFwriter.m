function [] = TIFwriter(img, savepath, compression)
    if nargin < 3
        img = uint16(img);
    
        % 2D image
        if length(size(img)) == 2
            imwrite(img, savepath, 'WriteMode', 'overwrite');
        end
    
        % 3D image
        if length(size(img)) == 3
            numberOfPages = size(img, 3);
            for frame_iter = 1:numberOfPages
                if frame_iter == 1
                    imwrite(img(:, :, frame_iter), savepath, 'WriteMode', 'overwrite');
                else
                    imwrite(img(:, :, frame_iter), savepath, 'WriteMode', 'append');
                end
            end    
        end
    else
        % 2D image
        if length(size(img)) == 2
            imwrite(img, savepath, 'WriteMode', 'overwrite', 'Compression', compression);
        end
    
        % 3D image
        if length(size(img)) == 3
            numberOfPages = size(img, 3);
            for frame_iter = 1:numberOfPages
                if frame_iter == 1
                    imwrite(img(:, :, frame_iter), savepath, 'WriteMode', 'overwrite', 'Compression', compression);
                else
                    imwrite(img(:, :, frame_iter), savepath, 'WriteMode', 'append', 'Compression', compression);
                end
            end    
        end
    end
end