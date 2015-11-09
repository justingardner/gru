% plot mystr2num
%
% usage: mystr2num
function ret = mystr2num(s)

if (isempty(s))
  ret = 0;
else
  ret = str2num(s);
  if (isempty(ret)) ret = 0; end
end


