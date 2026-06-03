function [row1, row2, col1, col2] = getROIboundary(roi_bw, roi_width)

    [rows, cols] = find(roi_bw);
    row_range = [min(rows), max(rows)]; col_range = [min(cols), max(cols)];
    h = size(roi_bw, 1); w = size(roi_bw, 2); 

    row1 = 1; row2 = roi_width; 
    col1 = 1; col2 = roi_width;

    if row_range(2)-row_range(1)+1 < roi_width
        if row_range(1)==1
            row1 = row_range(1)-row_range(2)+roi_width;
        elseif row_range(2)==h
            row2 = row_range(2)-row_range(1)+1;
        end
    end

    if  col_range(2)-col_range(1)+1 < roi_width
        if col_range(1)==1
            col1 = col_range(1)-col_range(2)+roi_width;
        elseif col_range(2)==w
            col2 = col_range(2)-col_range(1)+1;
        end
    end
end