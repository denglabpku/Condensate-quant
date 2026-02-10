function [img] = bfopenframe_oneChannel(image_file,frame_id)
%bfopenframe Open one frame from ND2 image with one channel
%   Detailed explanation goes here
stitchFiles = 0;
r = bfGetReader(image_file, stitchFiles);
numSeries = r.getSeriesCount();
r.setSeries(numSeries - 1);
img = bfGetPlane(r, frame_id);
end