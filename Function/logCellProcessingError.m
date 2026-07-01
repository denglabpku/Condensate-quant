function logCellProcessingError(error_log_path, filename, ME)
% Append the failed cell name and MATLAB error details to a text log.
fid = fopen(error_log_path, 'a');
if fid == -1
    warning('CellProcessing:LogOpenFailed', ...
        'Could not write failed cell to error log: %s', error_log_path);
    return
end

fprintf(fid, '[%s] %s\n', char(datetime('now')), filename);
fprintf(fid, 'Error ID: %s\n', ME.identifier);
fprintf(fid, 'Message: %s\n', ME.message);
for stack_iter = 1:numel(ME.stack)
    fprintf(fid, '  at %s line %d\n', ME.stack(stack_iter).name, ME.stack(stack_iter).line);
end
fprintf(fid, '\n');
fclose(fid);
end