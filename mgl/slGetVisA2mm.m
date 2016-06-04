

%slGetVisA2mm.m
%
%
% author: steeve laquitaine
%   date: 151203
%purpose: calculate visual angle distance (deg) in mm
%
%  usage:
%
%       dOnScreen = slGetVisA2mm(57,1)
%
% d2screen: eye distance to screen in mm
%     visA: visual angle in deg
%
%
%Description:
%- trigonometry
%
%Quick test: at 570 mm from the screen a distance of 1 deg
%            on screen should correspond to 1 mm.
%
%            dOnScreen = slGetVisA2mm(57,1)
%
%
%related: slGetScreenMMbyPix.m

function dOnScreen = slGetVisA2mm(d2screen,visA)

%convert deg to rad (unsigned ((1:2*pi)))
visArads = (visA/360)*2*pi;

%distance on screen in mm
dOnScreen = tan(visArads)*d2screen;
