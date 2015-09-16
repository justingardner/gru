% applyEventRelatedFiltering.m
%
%        $Id: applyEventRelatedFiltering.m,v 1.2 2008/07/16 21:38:18 justin Exp $ 
%      usage: timecourse = applyEventRelatedFiltering(timecourse,concatInfo)
%         by: justin gardner
%       date: 07/15/08
%    purpose: apply the same filtering that was done with the eventRelated
%             code to the timecourse. Timecourse is assumed to be a kxn
%             array where k is the number of time courses and n is the
%             number of time points
%
function timecourse = applyEventRelatedFiltering(timecourse,concatInfo)

% check arguments
if ~any(nargin == [1 2])
  help applyEventRelatedFiltering
  return
end

% if there is no concatInfo, just detrend and remove the mean
if ieNotDefined('concatInfo')
  % detrend
  timecourse = eventRelatedDetrend(timecourse')';
  % remove mean
  timecourse = timecourse-repmat(mean(timecourse')',1,size(timecourse,2));
  % add 1 back (so that it can be survive being converted back into % signal)
  timecourse = timecourse+1;
  return
else
  for runnum = 1:concatInfo.n
    % get the start and end vol for this part of the concatenation
    startVol = concatInfo.runTransition(runnum,1);
    endVol = concatInfo.runTransition(runnum,2);
    % detrend
    timecourse(startVol:endVol) = eventRelatedDetrend(timecourse(startVol:endVol));
    % use the hipass-filtering
    if isfield(concatInfo,'hipassfilter') && (length(concatInfo.hipassfilter) >= runnum)
      timecourse(startVol:endVol) = ifft(fft(timecourse(startVol:endVol)) .* concatInfo.hipassfilter{runnum});
    end
    % deal with projection
    if isfield(concatInfo,'projection') && ~isempty(concatInfo.projection{runnum})
      projectionWeight = concatInfo.projection{runnum}.sourceMeanVector * timecourse(startVol:endVol);
      timecourse(startVol:endVol) = timecourse(startVol:endVol) - concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
    end
    % remove mean
    timecourse(startVol:endVol) = timecourse(startVol:endVol)-repmat(mean(timecourse(startVol:endVol)')',1,endVol-startVol+1);
    % add 1 back (so that it can be survive being converted back into % signal)
    timecourse(startVol:endVol) = timecourse(startVol:endVol)+1;
  end
end
