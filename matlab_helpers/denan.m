% denan.m
%
%      usage: denan.m()
%         by: justin gardner
%       date: 04/10/03
%    purpose: removes nans from an array
function d = denan(d)

d = d(find(~isnan(d)));