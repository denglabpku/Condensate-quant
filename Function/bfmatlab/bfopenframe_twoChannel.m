function [Ch1_img,Ch2_img] = bfopenframe_twoChannel(image_file,frame_id)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
stitchFiles = 0;
r = bfGetReader(image_file, stitchFiles);
numSeries = r.getSeriesCount();
r.setSeries(numSeries - 1);
Ch1_img = bfGetPlane(r,2*(frame_id-1)+1);
Ch2_img = bfGetPlane(r,2*(frame_id-1)+2);
end