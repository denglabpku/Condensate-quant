function [CDboundary, labels] = getCDboundary(img, img_class, CDcenter, CDinterface, sep_point)
    %%%%%%%%%%%%%%%%%%%%%%%%%% CD boundary %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mask = img_class > sep_point;
    temp_labels = -ones(size(mask));

    temp = CDcenter & mask;
    ccomp = bwconncomp(temp);
    for comp_iter = 1:length(ccomp.PixelIdxList)
        temp_labels(ccomp.PixelIdxList{1, comp_iter}) = comp_iter;
    end
    
    % first round
    temp_img = img;
    temp_img(1, :)=0; temp_img(end, :)=0; temp_img(:, 1)=0; temp_img(:, end)=0; % remove boundary
    temp_img(~mask) = 0; temp_labels(temp_img==0) = 0;
    intensity = temp_img(:);
    [~, sort_idx] = sort(intensity, 'descend');
    for i = 1:length(sort_idx)
        idx = sort_idx(i);
        if intensity(idx) > 0
            yid = mod((idx - 1), size(temp_img, 1)) + 1;
            xid = (idx - yid)/size(temp_img, 1) + 1;
            neighbor_x = [xid-1; xid; xid+1; xid-1; xid+1; xid-1; xid; xid+1];
            neighbor_y = [yid-1; yid-1; yid-1; yid; yid; yid+1; yid+1; yid+1]; % 8-connected
            neighbor_idx = neighbor_y + (neighbor_x-1)*size(temp_img, 1);
            if temp_labels(idx) == -1
                if sum(temp_labels(neighbor_idx) > 0) > 0
                    temp = [temp_labels(neighbor_idx), intensity(neighbor_idx)];
                    temp = temp(temp(:, 1)>0, :);
                    [~, temp_idx] = max(temp(:, 2));
                    temp_labels(idx) = temp(temp_idx, 1);
                end
            end
        else
            break
        end
    end
    
    % second round
    temp_img = img;
    temp_img(temp_labels~=-1) = 0;
    intensity = temp_img(:);
    [~, sort_idx] = sort(intensity, 'descend');
    currentIdx = length(ccomp.PixelIdxList);
    for i = 1:length(sort_idx)
        idx = sort_idx(i);
        if intensity(idx) > 0
            yid = mod((idx - 1), size(temp_img, 1)) + 1;
            xid = (idx - yid)/size(temp_img, 1) + 1;
            neighbor_x = [xid-1; xid; xid+1; xid-1; xid+1; xid-1; xid; xid+1];
            neighbor_y = [yid-1; yid-1; yid-1; yid; yid; yid+1; yid+1; yid+1]; % 8-connected
            neighbor_idx = neighbor_y + (neighbor_x-1)*size(temp_img, 1);
            if temp_labels(idx) == -1
                if sum(temp_labels(neighbor_idx) > 0) > 0
                    aa = temp_labels(neighbor_idx); aa = aa(aa>0);
                    temp_labels(idx) = mode(aa);
                else
                    currentIdx = currentIdx + 1;
                    temp_labels(idx) = currentIdx;
                end
            end
        else
            break
        end
    end
    
    % third round
    unassigned = temp_labels>length(ccomp.PixelIdxList);
    ccomp1 = bwconncomp(unassigned);
    for comp_iter = 1:length(ccomp1.PixelIdxList)
        temp_ccomp = logical(zeros(size(temp_labels)));
        temp_ccomp(ccomp1.PixelIdxList{1, comp_iter}) = 1;
        bb = temp_labels(logical(imdilate(temp_ccomp, strel('disk', 1)) - temp_ccomp));
        bb = bb(bb>0);
        if ~isempty(bb)
            temp_labels(ccomp1.PixelIdxList{1, comp_iter}) = mode(bb);
        else
            temp_labels(ccomp1.PixelIdxList{1, comp_iter}) = 0;
        end
    end

    labels = temp_labels;

    temp_img_class = zeros(size(temp_img, 1)+2, size(temp_img, 2)+2);
    temp_img_class(2:(end-1), 2:(end-1)) = temp_labels;
    temp_CDboundary = zeros(size(temp_img));
    for x = 2:(size(temp_img_class, 2)-1)
        for y = 2:(size(temp_img_class, 1)-1)
            if temp_img_class(y, x) > 0
                temp_neighbor = temp_img_class((y-1):(y+1), (x-1):(x+1));
                if sum(temp_neighbor(:) > temp_img_class(y, x)) > 0
                    temp_CDboundary(y-1, x-1) = 1;
                end
            end
        end
    end
    CDboundary = logical(temp_CDboundary) | CDinterface;

end

