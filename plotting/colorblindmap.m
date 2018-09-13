function [cmap,names] = colorblindmap()
%% COLORBLINDMAP
% A set of colors that are unambiguous to common color blind users. The
% first index is black. Colors are in RGB.
%
% Dan Birman 2018/05/30
%
% From http://jfly.iam.u-tokyo.ac.jp/color/

cmap = [0,0,0
          230,159,0
          86,180,233
          0,158,115
          240,228,66
          0,114,178
          213,94,0
          204,121,167];
names = {'Black','Orange','Sky Blue','Bluish Green','Yellow','Blue','Vermillion','Pink'};