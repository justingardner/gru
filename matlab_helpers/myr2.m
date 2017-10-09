function r2 = myr2( y,y_ )
% R^2
%
% Use:
%   r2 = myr2(y,y')
% Input: 
%   y  - empirical data
%   y' - predicted data
%
% Output:
%   r2 - R^2 calculated according to: 1 - SSE/SST
%
% WARNING: uses correlation coefficent ^2

if ~iscolumn(y)
    y = y(:);
end

if ~iscolumn(y_)
    y_ = y_(:);
end

% r2 = 1 - (sum((y_-y).^2)/sum((y-mean(y)).^2));
cc = corrcoef([y y_]);
r2 = cc(1,2)^2;

