% weibull.m
%
%      usage: weibullInv(x,params,<logBase>)
%         by: justin gardner
%       date: 07/29/10
%    purpose: inverse of weibull function, params are [T beta lambda gamma]
%             T is the threshold
%             beta is the slope
%             lambda is the maximum value of the function (i.e. 
%             the lapse rate - or the rate at which the observer
%             lapses or makes "finger errors")
%             gamma is the minimum value of the function (i.e. when
%             the signal strength is low and the subject is
%             guessing-should be 0.5 for 2AFC)
%
%             You can also compute the weibull where the input values (x) are
%             log values. For instance, if x are log base 10 units, then
% 
%             x = weibullInv(y,params,10);
% 
%             Note that then, threshold parameter should be specified in log base 10. (beta, gamma and delta are the same)
% 
function retval = weibullInv(y,params,logBase)

if (nargin == 1)
  params = [2 3.5 0.01 0.5];
elseif (nargin == 2)
  logBase = [];
elseif nargin ~= 3
  help weibull
  return
end

% threshold
T = params(1);
% slope of psychometric function
beta = params(2);
% lapse rate
lambda = params(3);
% probability of response at 0
gamma = params(4);

y=(y-gamma)/(1-lambda-gamma);

% the weibull function
if isempty(logBase)
  retval = wblinv(y,T,beta);
else
  retval = wblinv(y,logBase.^T,beta);
end  




