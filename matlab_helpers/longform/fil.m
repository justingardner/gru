
function data = fil(data,col,comparison,value)
% LONGFORM FILTER(data,column,>/</>=/<=,==,!=,value)

scol = data(:,col);
switch comparison
    case '>'
        data = data(scol>value,:);
    case '<'
        data = data(scol<value,:);
    case '>='
        data = data(scol>=value,:);
    case '<='
        data = data(scol<=value,:);
    case '=='
        data = data(scol==value,:);
    case '!='
        data = data(scol~=value,:);
end