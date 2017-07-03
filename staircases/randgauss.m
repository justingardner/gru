% randgauss.m
%
%      usage: randgauss(mean,std,n)
%         by: justin gardner
%       date: 04/13/04
%    purpose: create random gaussian distribution
%       e.g.: 
%
function randdist = randgauss(mean,std,n)

% check input arguments
if (nargin ~= 3)
  help randgauss;
  return
end

% make a uniform distribution
randnums = rand(1,n);
% pass through inverse cumulative gaussian function
% and mutiply by desired std and add mean
randdist = sqrt(2)*std*erfinv(randnums*2-1)+mean;
return

% old way of doing this (required for loop...sloooooww.

parameters = [0 1/sqrt(2) 1];
randnums = rand(1,n);
for i = 1:n
    randdist2(i) = icg(parameters,randnums(i));
end

%%%%%%%%%%%%%%%%%%%%%%%
% cumulative gaussian
%%%%%%%%%%%%%%%%%%%%%%%
function retval = cg(parameters, xpoints)

% decode parameters
mean = parameters(1);
sd = parameters(2);
max = parameters(3);

% calculate cumulative gaussian
% erf goes from -1 to 1 so we have to scale
% properly. Also, we are scaling the input points
% so we don't have to scale the output of the erf by
% sqrt(2)*sd. 
retval = (max)*((erf((xpoints-mean)/(sqrt(2)*sd))+1)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse cumulative gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = icg(parameters, ypoint)

min = -500; max = 500;
epsilon = .00001;

% return immediately with infinity if ypoint is
% beyond range of curve
if (parameters(3) < ypoint)
  retval = inf;
  return;
end

err = inf;
% search by bisection for x that matches ypoint
while (abs(err) > epsilon)
  % bisect interval
  mid = (max-min)/2+min;

  % figure out which side the ypoint lives in
  midy = cg(parameters,mid);

  % calculate error between point and desired ypoint
  err = ypoint - midy;

  % move interval accordingly
  if (ypoint > midy)
    min = mid;
  else
    max = mid;
  end
  %disp(sprintf('%0.5f:%0.5f',min,max));
end
retval = mid;

