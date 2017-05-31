% motionEnergyModel.m
%
%        $Id:$ 
%      usage: r = motionEnergyModel(s)
%         by: justin gardner
%       date: 05/29/17
%    purpose: implements motion-energy (Adelson & Bergen, 1985) model 
%             and compute responses on passed in stimulus
%
%             Read code comments for settable parameters of model
%
%             [s msc] = motionEnergyModelMakeStimulus('offscreen');
%             r = motionEnergyModel(s,'myscreen',msc);
function r = motionEnergyModel(s,varargin)

% check arguments
if nargin < 1
  help motionEnergyModel
  return
end

% parse args
getArgs(varargin,{'stimulusSize=5','orientationPreference',0:45:359,'sfPreference',[1],'tfPreference',[1],'maxFilterTime',1,'deltaTime=0.1','deltaSpace=0.1','dispFigures',1,'myscreen',[]});

% get parameters from myscreen
if ~isempty(myscreen)
  [deltaTime stimulusSize deltaSpace] = getParametersFromMyscreen(myscreen);
end

% make displays for debugging etc
m.dispFigures = dispFigures;

% specify model (for now we just have adelson-bergen)
m.linearSpatialTemporalFilter = 'adelson-bergen';

% time steps for rf computation (in seconds)
m.deltaTime = deltaTime; 

% length of time to compute receptive field for (in seconds)
m.maxFilterTime = maxFilterTime;

% spatial steps for rf computation (in deg)
m.deltaSpace = deltaSpace;

% size of stimulus in degrees which will determine the
% max size of all filters - i.e. the filters will all be
% calculated to fill out the whole stimulus
% The linear receptive field within this area 
% is always centered within this and that the "size" of the 
% rf will largely be determined by the desired spatial frequency.
% May want to update this code to have rfs that fit 
% completely within the window to save processing time
m.stimulusSize = stimulusSize;

% orientations of filter in degrees
m.orientationPreference = orientationPreference;

% the spatial frequencies over which to compute filters in cycles/deg
m.sfPreference = sfPreference;

% the temporal frequencies over which to compute filters in cycles/sec
m.tfPreference = tfPreference;

for iTF = 1:length(m.tfPreference)
  for iSF = 1:length(m.sfPreference)
    % compute necessary sigma and k for this sf/tf
    [m.sigma(iSF) m.k(iTF)] = getSigmaK(m,m.sfPreference(iSF),m.tfPreference(iTF));

    % compute linear receptive
    m.baseFilter = computeSpatialTemporalFilter(m,m.sigma(iSF),m.k(iTF));

    % rotate to get different orientations
    m.filters(iTF,iSF,:) = getDifferentOrientations(m,m.baseFilter(1),m.orientationPreference);
  end
end

% normalize all filters to have the same total sum-of-squares
m.filters = normalizeLinearFilters(m.filters);

% apply the filters to the stimulsu
r = applyFilters(s,m);

% display what we got
dispFilterResponse(m,r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispFilterResponse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispFilterResponse(m,r)

for iTF = 1:length(m.tfPreference)
  for iSF = 1:length(m.sfPreference)
    for iOrient = 1:length(orientationPreference)
      
    end
  end
end
keyboard

%%%%%%%%%%%%%%%%%%%%%%
%    applyFilters    %
%%%%%%%%%%%%%%%%%%%%%%
function r = applyFilters(s,m)

disppercent(-inf,'(motionEnergyModel) Applying filters to stimulus');
for iFilter = 1:length(m.filters(:))
  phase1 = myTimeConv(s,m.filters(iFilter).phase1);
  phase2 = myTimeConv(s,m.filters(iFilter).phase2);
  r(iFilter,:) = sqrt(phase1.^2 + phase2.^2);
  disppercent(iFilter/length(m.filters(:)));
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%    myTimeConv    %
%%%%%%%%%%%%%%%%%%%%
function r = myTimeConv(s,f)

% get length of stimulus and filter
slen = size(s,3);
flen = size(f,3);

% preallocate r
r = nan(size(s,1),size(s,2),max([slen+flen-1,slen,flen]));

% do 1D convolution at each location in space (surely there is
% a faster way to do this?
for iX = 1:size(s,1)
  for iY = 1:size(s,2)
    r(iX,iY,:) = conv(squeeze(s(iX,iY,:)),squeeze(f(iX,iY,:)));
  end
end

% sum over space
r = squeeze(sum(sum(r)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    normalizeLInearFilters    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filters = normalizeLinearFilters(filters)

% loop over filters
for iFilter = 1:length(filters(:))
  % normalize by sum of squares
  filters(iFilter).phase1 = filters(iFilter).phase1/sqrt(sum(filters(iFilter).phase1(:).^2));
  filters(iFilter).phase2 = filters(iFilter).phase2/sqrt(sum(filters(iFilter).phase2(:).^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getDifferentOrientations    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filters = getDifferentOrientations(m,baseFilter,orientationPreferences)

% get x and y meshgrids
x = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;
y = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;
[xGrid yGrid] = meshgrid(x,y);

% number of orientations
nOrient = length(orientationPreferences);

disppercent(-inf,'(motionEnergyModel) Rotating to get different orientation/direction selectivity');
for iOrient = 1:nOrient
  % make rotation matrix
  c = cos(d2r(orientationPreferences(iOrient)));
  s = sin(d2r(orientationPreferences(iOrient)));
  rotationMatrix = [c -s ; s c];
  % rotate coordinates
  rotatedCoordinates = [xGrid(:) yGrid(:)] * rotationMatrix;
  xRotatedGrid = reshape(rotatedCoordinates(:,1),length(x),length(y));
  yRotatedGrid = reshape(rotatedCoordinates(:,2),length(x),length(y));
  % interpolate across time
  for iTime = 1:size(baseFilter.phase1,3)
    % rotate phase 1
    filters(iOrient).phase1(:,:,iTime) = interp2(xGrid,yGrid,squeeze(baseFilter.phase1(:,:,iTime)),xRotatedGrid,yRotatedGrid,'linear',0);
    % rotate phase 2
    filters(iOrient).phase2(:,:,iTime) = interp2(xGrid,yGrid,squeeze(baseFilter.phase2(:,:,iTime)),xRotatedGrid,yRotatedGrid,'linear',0);
  end
  disppercent(iOrient/nOrient);
end
disppercent(inf);

% now display 
if m.dispFigures
  % init figure
  mlrSmartfig('moitonEnergyModel_orientations');
  
  % compute maxFilterTimePoint
  for iTime = 1:size(baseFilter.phase1,3)
    timeEnergy(iTime) = sqrt(sum(sum(baseFilter.phase1(:,:,iTime).^2)));
  end
  [~,maxFilterTimePoint] = max(timeEnergy);
  for iOrient = 1:nOrient
    subplot(nOrient,3,((iOrient-1)*3)+1);
    dispXY(filters(iOrient).phase1,maxFilterTimePoint);axis off
    subplot(nOrient,3,((iOrient-1)*3)+2);
    dispXT(filters(iOrient).phase1);axis off
    subplot(nOrient,3,((iOrient-1)*3)+3);
    dispYT(filters(iOrient).phase1);axis off
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeSpatialTemporalFilter    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filter = computeSpatialTemporalFilter(m,sigma,k)

filter = [];

% this function is just a switch for different filter types
switch(fixBadChars(lower(m.linearSpatialTemporalFilter))) 
 % adelson bergen
 case {fixBadChars(lower('adelson-bergen'))}
   filter = computeSpatialTemporalFilterAdelsonBergen(m,sigma,k);
end

%%%%%%%%%%%%%%%%%%%
%    getSigmaK    %
%%%%%%%%%%%%%%%%%%%
function [sigma k] = getSigmaK(m,sf,tf)

% make sure sample time is sufficient for this
m.deltaTime = 0.001;
m.maxFilterTime = 10;

% get best fit k
k = fminsearch(@tfMin,100,optimset,m,tf);

% check if we have significant error in the best TF
if tfMin(k,m,tf)>0.01
  disp(sprintf('(motionEnergyModel) Could not achieve peak filter response at %f',tf));
  keyboard
end

% make sure sample space is sufficient for this
m.deltaSpace = 0.0001;
m.stimulusSize = 10;

% get best fit sigma
sigma = fminsearch(@sfMin,0.01,optimset,m,sf);

% check if we have significant error in the best TF
if sfMin(sigma,m,sf)>0.01
  disp(sprintf('(motionEnergyModel) Could not achieve peak filter response at %f',sf));
  keyboard
end


if m.dispFigures
  % compute temporal frequency response and peak TF (see tfMin)
  t = 0:m.deltaTime:m.maxFilterTime;
  temporalResponse = temporalImpulseResponse(t,k,3);
  temporalResponseFFT = abs(fft(temporalResponse));
  [~,maxTFBin] = max(temporalResponseFFT);
  TFpeak = (maxTFBin-1)/(m.maxFilterTime+m.deltaTime);
  
  % compute spatial frequency response and peak TF (see sfMin)
  x = [0 1];
  y = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;
  [xGrid,yGrid] = meshgrid(x,y);
  [sf1 sf2] = spatialProfile(xGrid,yGrid,sigma);
  sf1FFT = abs(fft(sf1(:,1)));
  [~,maxSFBin] = max(sf1FFT);
  SFpeak = (maxSFBin-1)/(m.stimulusSize+m.deltaSpace);

  % plot temporal frequency response / fft
  mlrSmartfig('getSigmaK','reuse');clf;
  subplot(2,2,1);
  plot(t,temporalResponse);
  subplot(2,2,2);
  n = maxTFBin*4;
  plot((0:n-1)/(m.maxFilterTime+m.deltaTime),temporalResponseFFT(1:n));
  xlabel('Temporal frequency (cycles / sec)');
  vline(TFpeak);
  title(sprintf('FFT: max = %0.2f c/sec',TFpeak));

  % plot spatial frequency response / fft
  subplot(2,2,3);
  plot(y,sf1(:,1));
  subplot(2,2,4);
  n = maxSFBin*4;
  plot((0:n-1)/(m.stimulusSize+m.deltaSpace),sf1FFT(1:n));
  vline(SFpeak);
  title(sprintf('FFT: max = %0.2f c/deg',SFpeak));
  drawnow
end

%%%%%%%%%%%%%%%
%    tfMin    %
%%%%%%%%%%%%%%%
function r = tfMin(k,m,tf)

t = 0:m.deltaTime:m.maxFilterTime;

% get max of fft of impulse response
t1 = temporalImpulseResponse(t,k,3);
t1 = t1-mean(t1);
temporalResponseFFT = abs(fft(t1));
[~,maxBin] = max(temporalResponseFFT(2:end));

% convert to tf
TFpeak = maxBin/(m.maxFilterTime+m.deltaTime);

% get difference from desired
r = (TFpeak-tf).^2;

%%%%%%%%%%%%%%%
%    sfMin    %
%%%%%%%%%%%%%%%
function r = sfMin(sigma,m,sf)

% get spatial profile (note that we have
% to have x have more than one value so that
% we can use the gradient function in spatialProfile)
x = [0 1];
y = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;

[xGrid,yGrid] = meshgrid(x,y);
[sf1 sf2] = spatialProfile(xGrid,yGrid,sigma);
sf1 = sf1(:,1);
sf1 = sf1-mean(sf1);

% get max of fft of sf1
sf1fft = abs(fft(sf1));
[~,maxBin] = max(sf1fft(2:end));

% convert to sf
SFpeak = maxBin/(m.stimulusSize+m.deltaSpace);

% get difference from desired
r = (SFpeak-sf).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    temporalImpulseResponse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = temporalImpulseResponse(t,k,n)

r = ((k*t).^n) .* exp(-k*t) .* (1/factorial(n) - ((k*t).^2)./factorial(n+2));

%%%%%%%%%%%%%%%%%%%%%%%%
%    spatialProfile    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [sf1 sf2] = spatialProfile(x,y,sigma)

% compute gaussian receptive field
sf = exp (-(x.^2 + y.^2) / (sigma^2 + sigma^2));

% sf1 is the 2nd derivative
[~,sf1] = gradient(sf);
[~,sf1] = gradient(sf1);
sf1 = -sf1;
% normalize
sf1 = sf1/sqrt(sum(sf1(:).^2));

% sf2 is the 3rd derivative
[~,sf2] = gradient(sf1);
sf2 = sf2/sqrt(sum(sf2(:).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeSpatialTemporalFilterAdelsonBergen    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filter = computeSpatialTemporalFilterAdelsonBergen(m,sigma,k)

% amount of time to compute response for
t = 0:m.deltaTime:m.maxFilterTime;

% compute the TF response of Adelson & Bergen (Equation 1) - the two TF functions
% have n = 3 and n = 5 as in text.
tf1 = temporalImpulseResponse(t,k,3);
tf2 = temporalImpulseResponse(t,k,5);

% compute the SF response
% first, the x and y dimensions
x = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;
y = -m.stimulusSize/2:m.deltaSpace:m.stimulusSize/2;
[xGrid yGrid] = meshgrid(x,y);

% now compute the spatial receptive field
[sf1 sf2] = spatialProfile(xGrid,yGrid,sigma);

% now make the four combinations of sf and tf filters
rf11 = zeros(length(x),length(y),length(t));
rf12 = zeros(length(x),length(y),length(t));
rf21 = zeros(length(x),length(y),length(t));
rf22 = zeros(length(x),length(y),length(t));

disppercent(-inf,'(motionEnergyModel) Making combinations of SF and TF filters');
for iX = 1:length(x)
  for iY = 1:length(y)
    rf11(iX,iY,:) = sf1(iX,iY)*tf1;
    rf12(iX,iY,:) = sf1(iX,iY)*tf2;
    rf21(iX,iY,:) = sf2(iX,iY)*tf1;
    rf22(iX,iY,:) = sf2(iX,iY)*tf2;
  end
  disppercent(iX/length(x));
end
disppercent(inf);

% now compute the sums and differences.
filter(1).dir = 0;
filter(1).phase1 = rf21 + rf12;
filter(1).phase2 =  rf11 - rf22; 

filter(2).dir = 180;
filter(2).phase1 = rf12 - rf21;
filter(2).phase2 = rf22 + rf11; 


% make plot
if m.dispFigures
  % init figure
  mlrSmartfig('Adelson_Bergen_1','reuse');clf;
  % get middle point of spatial rf
  midPoint = round(size(sf1,1)/2);
  % display slice of spatial RF for both filters
  subplot(4,4,1:4);
  plot(x,sf1(:,midPoint));
  hold on
  plot(x,sf2(:,midPoint));
  hline(0);
  ylabel('Filter Response');
  xlabel('Space (deg)');
  title('Spatial slice of impluse response');
  % display temporal RF for both filters
  subplot(4,4,5:8);
  plot(t,tf1);
  hold on
  plot(t,tf2);
  hline(0);
  xlabel('Time (sec)');
  ylabel('Filter Response');
  title('Temporal response');
  % display slices of images
  subplot(4,4,9)
  dispXT(rf21);
  subplot(4,4,10)
  dispXT(rf22)
  subplot(4,4,11)
  dispXT(rf12);
  subplot(4,4,12)
  dispXT(rf11)

  subplot(4,4,13)
  dispXT(filter(1).phase1);
  subplot(4,4,14)
  dispXT(filter(1).phase2);
  subplot(4,4,15)
  dispXT(filter(2).phase1);
  subplot(4,4,16)
  dispXT(filter(2).phase2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    display image with correct orientation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispXT(im)

midpoint = round(size(im,2)/1);
imagesc(fliplr(squeeze(im(:,midpoint,:))'));
colormap(gray);
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    display image with correct orientation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispYT(im)

midpoint = round(size(im,1)/1);
imagesc(fliplr(squeeze(im(midpoint,:,:))'));
colormap(gray);
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    display image with correct orientation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispXY(im,t)

imagesc(fliplr(squeeze(im(:,:,t))'));
colormap(gray);
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getParametersFromMyscreen    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deltaTime stimulusSize deltaSpace] = getParametersFromMyscreen(myscreen)


% check that the image is square
if ~isequal(myscreen.imageWidth,myscreen.imageHeight)
  disp(sprintf('(motionEnergyModel) No support yet for non-square screens'));
  keyboard
end

% get parameters
deltaTime = myscreen.frametime;
stimulusSize = myscreen.imageWidth;
deltaSpace = stimulusSize/(myscreen.screenWidth-1);

