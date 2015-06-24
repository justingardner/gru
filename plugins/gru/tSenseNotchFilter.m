% tSenseNotchFilter
%
%      usage: [tSeries notchFilter] = tSenseNotchFiler(tSeries,params)
%         by: justin gardner
%       date: 11/21/2012
%    purpose: Runs tSense notch filter on data. The last dimension of data is the
%             time series (so, it can be an 1xn, ixjxkxn or whatever).
%             params can be either the tSense acceleration factor or a structure
%             (when called from tSenseNotch that has the filed tSenseAcc).
%
function [tSeries notchFilter] = tSenseNotchFilter(tSeries,params)

notchFilter = [];
if ~any(nargin == [2])
  help tSenseNotchFilter
  tSeries = [];
  return
end

% interp params if it is set as a number
if isnumeric(params)
  acc = params;
  params = [];
  params.tSenseAcc = acc;
end

% set default filterHalfWidth (0 mean knock out single frequency)
if ~isfield(params,'notchFilterHalfWidth')
  params.notchFilterHalfWidth = 0;
end
if ~isfield(params,'dispFilter')
  params.dispFilter = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the tSense notch filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first get length of filter
dataDim = size(tSeries);
n = dataDim(ndims(tSeries));
notchFilter = ones(1,n);


if iseven(n)
  if params.tSenseAcc == 2
    notchFilter((n/2)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    disp(sprintf('(tSenseNotchFilter) Applying notch filter for acceleration of x%i',params.tSenseAcc));
  elseif params.tSenseAcc == 4
    disp(sprintf('(tSenseNotchFilter) Using filter width of %i',params.notchFilterHalfWidth));
    notchFilter((n/2)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    notchFilter((n/4)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    notchFilter((3*n/4)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    disp(sprintf('(tSenseNotchFilter) Applying notch filter for acceleration of x%i',params.tSenseAcc));
  elseif params.tSenseAcc == 3
    notchFilter(floor(n/3)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    notchFilter((2*floor(n/3))+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    disp(sprintf('(tSenseNotchFilter) Applying notch filter for acceleration of x%i',params.tSenseAcc));
  else
    mrErrorDlg(sprintf('(tSenseNotchFilter) Notch for tSense not implemented yet for acceleration of %i',params.tSenseAcc));
  end
else
  if params.tSenseAcc == 3
    notchFilter(floor(n/3)+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    notchFilter((2*floor(n/3))+1+(-params.notchFilterHalfWidth:params.notchFilterHalfWidth)) = 0;
    disp(sprintf('(tSenseNotchFilter) Applying notch filter for acceleration of x%i',params.tSenseAcc));
  else
    mrErrorDlg(sprintf('(tSenseNotchFilter) Notch for tSense not implemented yet for acceleration of %i',params.tSenseAcc));
  end
end

% display the filter if called for
if (params.dispFilter)
  mlrSmartfig('tSenseNotchFilter','reuse');clf
  subplot(1,2,1);
  plot(notchFilter,'k.-');
  hold on
  zoom on
  xlabel('frequency component');
  ylabel('magnitude');
  xaxis(1,n);
  yaxis(-0.2,1.2);
  title(sprintf('Notch filter for: Acceleration %i with half-width %i (0 means single component notch)',params.tSenseAcc,params.notchFilterHalfWidth));

  
  % display the convolution kernel of the filter
  subplot(1,2,2);
  plot(abs(fft(notchFilter)));
  xlabel('time (units uncalculated)');
  ylabel('filter magnitude');
  xaxis(0,n);
end

% faster way, do across first dimension
% replicate the filter
for i = 1:dataDim(1)
  newNotchFilter(i,:) = notchFilter;
end
% filter
for k = 1:dataDim(3)
  disppercent(k/dataDim(3));
  for j = 1:dataDim(2)
    timecourses = squeeze(tSeries(:,j,k,:));
    timecourses = ifft(fft(timecourses,[],2) .* newNotchFilter,[],2);
    tSeries(:,j,k,:) = real(timecourses);
  end
end
disppercent(inf);



 
