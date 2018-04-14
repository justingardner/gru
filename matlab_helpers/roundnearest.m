function vals = roundnearest( vals, nearvals )
%ROUNDNEAREST Round values to the nearest specified value
%
% Returns the first argument rounded to the nearest of the values in the
% second argument

for vi = 1:length(vals)
    d = abs(vals(vi)-nearvals);
    [~,di] = min(d);
    vals(vi) = nearvals(di);
end

