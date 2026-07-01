function initializeCellProcessingErrorLog(error_log_path)
% Start a fresh per-run cell processing error log.
fid = fopen(error_log_path, 'w');
if fid == -1
    warning('CellProcessing:LogOpenFailed', ...
        'Could not initialize error log: %s', error_log_path);
    return
end

fprintf(fid, 'Failed cell processing log\n');
fprintf(fid, 'Started: %s\n\n', char(datetime('now')));
fclose(fid);
end