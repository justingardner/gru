
%slGetScreenMMbyPix.m
%
%author: steeve laquitaine
%  date: 151203
%purpose: get mm by pixel on screen
%
%usage:
%
%       mmByPix = slGetScreenMMbyPix
%
%require mgl (justin gardner) library


function mmByPix = slGetScreenMMbyPix

%convert distance (in pix) to visual angle (deg)
dinfo = mglDescribeDisplays;
%screen true width in mm
widthmm = dinfo.screenSizeMM(1);
%screen true height in mm
heightmm = dinfo.screenSizeMM(2);
%screen width in pix
widthpix = dinfo.screenSizePixel(1);
%screen height in pix
heightpix = dinfo.screenSizePixel(2);
%mm per pix
mmByPix = min([widthmm/widthpix heightmm/heightpix]);
