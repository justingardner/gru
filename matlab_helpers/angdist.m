
function d = angdist(t1,t2)
% compute the angle between two values
%
% usage
%   angle distance = angdist(angle1, angle2)

d = acos(cos(t1).*cos(t2)+sin(t1).*sin(t2));