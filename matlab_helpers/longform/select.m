
function data = select(data,col,value)

scol = data(:,col);
data = data(scol==value,:);