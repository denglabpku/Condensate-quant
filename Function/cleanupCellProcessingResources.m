function cleanupCellProcessingResources(reader)
% Close resources that may be left open when a cell fails midway.
if nargin > 0 && ~isempty(reader)
    try
        reader.close();
    catch
    end
end

try
    close(findall(0, 'Type', 'figure', 'Tag', 'TMWWaitbar'));
catch
end

try
    close all;
catch
end
end