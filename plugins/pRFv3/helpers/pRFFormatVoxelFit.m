function displayString = pRFFormatVoxelFit(fit,coords,displayPrefix,rfType,prefitOnly)
% pRFFormatVoxelFit
%
% Return the model-specific one-line voxel summary used by pRFFit verbose
% output and by the sampled pRF quality-control display. Keeping the format
% in one function prevents the two displays from drifting apart.

if nargin < 2 || isempty(coords),coords = [NaN NaN NaN];end
if nargin < 3 || isempty(displayPrefix),displayPrefix = '';end
if nargin < 4 || isempty(rfType)
  if isstruct(fit) && isfield(fit,'rfType')
    rfType = fit.rfType;
  else
    rfType = '';
  end
end
if nargin < 5 || isempty(prefitOnly)
  prefitOnly = isstruct(fit) && isfield(fit,'bestFitVoxel');
end

coords = double(coords(:)');
if numel(coords) < 3
  coords(end+1:3) = NaN;
end

% pRFFit historically does not print a result line when a fit is empty.
% QC samples still need an explicit indication that their selected fit did
% not return; successful fits below use the exact verbose representation.
if isempty(fit)
  displayString = sprintf('%s[%2.f %2.f %2.f] FIT FAILED/EMPTY', ...
    displayPrefix,coords(1),coords(2),coords(3));
  return
end

r2 = fitField(fit,'r2',NaN);
polarAngleDegrees = radiansToDegrees(fitField(fit,'polarAngle',NaN));
eccentricity = fitField(fit,'eccentricity',NaN);
rfHalfWidth = fitField(fit,'std',NaN);

% A prefit-only result returns before pRFFit's model-specific verbose
% switch, so it always uses the base Gaussian-style line.
if prefitOnly
  displayString = sprintf( ...
    '%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f', ...
    displayPrefix,coords(1),coords(2),coords(3),r2,polarAngleDegrees, ...
    eccentricity,rfHalfWidth);
  return
end

switch rfType
  case 'gaussian-css'
    displayString = sprintf( ...
      '%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f exp=%f', ...
      displayPrefix,coords(1),coords(2),coords(3),r2,polarAngleDegrees, ...
      eccentricity,rfHalfWidth,fitField(fit,'css',NaN));
  case 'divNorm'
    displayString = sprintf( ...
      ['%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f ' ...
       'stdCenter=%6.1f stdSurround=%6.1f gainCenter=%6.1f gainSurround=%6.1f actC=%6.1f divC=%6.1f'], ...
      displayPrefix,coords(1),coords(2),coords(3),r2,polarAngleDegrees, ...
      eccentricity,rfHalfWidth,fitField(fit,'stdSurround',NaN), ...
      fitField(fit,'gainCenter',NaN),fitField(fit,'gainSurround',NaN), ...
      fitField(fit,'actConst',NaN),fitField(fit,'divConst',NaN));
  otherwise
    displayString = sprintf( ...
      '%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f', ...
      displayPrefix,coords(1),coords(2),coords(3),r2,polarAngleDegrees, ...
      eccentricity,rfHalfWidth);
end


function value = fitField(fit,fieldName,defaultValue)

if isstruct(fit) && isfield(fit,fieldName)
  value = fit.(fieldName);
else
  value = defaultValue;
end


function degrees = radiansToDegrees(angle)

% Match pRFFit's historical local r2d helper exactly.
degrees = (angle/(2*pi))*360;
while sum(degrees > 360)
  degrees = degrees-(degrees > 360)*360;
end
while sum(degrees < -360)
  degrees = degrees+(degrees < -360)*360;
end
