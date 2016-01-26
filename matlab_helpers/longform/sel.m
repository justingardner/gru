
function data = sel(data,col,value)
% LONGFORM SELECT(data,column,==value)

scol = data(:,col);
data = data(scol==value,:);