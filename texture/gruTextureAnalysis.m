function [params intermediateSteps] = gruTextureAnalysis(im0, Nsc, Nor, Na)
%
% This function is taken from Portilla & Simoncelli code and just adds
% some comments to make it easier to read and stores some intermediate%
% steps for display by dispPSStats. Note that original comments are
% typically prefixed by %% and new comments are always on their own 
% line and are prefixed by %
%
% Analyze texture for application of Portilla-Simoncelli model/algorithm.
%
% [params] = textureAnalysis(im0, Nsc, Nor, Na);
% 	im0: 	original image
% 	Nsc: 	number of scales
% 	Nor: 	number of orientations
% 	Na:	spatial neighborhood considered (Na x Na)	
%
% Example: Nsc=4; Nor=4; Na=7;
%
% See also textureSynthesis.

% Javier Portilla and Eero Simoncelli.
% Work described in:
%  "A Parametric Texture Model based on Joint Statistics of Complex Wavelet Coefficients".
%  J Portilla and E P Simoncelli. Int'l Journal of Computer Vision,
%  vol.40(1), pp. 49-71, Dec 2000.   
%
% Please refer to this publication if you use the program for research or
% for technical applications. Thank you.
%
% Copyright, Center for Neural Science, New York University, January 2001.
% All rights reserved.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Warn = 0;  % Set to 1 if you want to see warning messages

%% Check required args are passed
if (nargin < 4)
  error('Function called with too few input arguments');
end

%% 1D interpolation filter, for scale cross-correlations:
interp = [-1/16 0 9/16 1 9/16 0 -1/16]/sqrt(2);

if ( mod(Na,2) == 0 )
  error('Na is not an odd integer');
end

%% If the spatial neighborhood Na is too big for the lower scales,
%% "modacor22.m" will make it as big as the spatial support at
%% each scale:

[Ny,Nx] = size(im0);
nth = log2(min(Ny,Nx)/Na);
if nth<Nsc & Warn,
  fprintf(1,'Warning: Na will be cut off for levels above #%d !\n', floor(nth+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

la = floor((Na-1)/2);

%% Pixel statistics
[mn0 mx0] = range2(im0);
mean0 = mean2(im0);
var0 = var2(im0, mean0);
skew0 = skew2(im0, mean0, var0);
kurt0 = kurt2(im0, mean0, var0);
statg0 = [mean0 var0 skew0 kurt0 mn0 mx0];

% Add a little bit of noise to the original, in case it has been 
% artificially generated, to avoid instability crated by symmetric
% conditions at the synthesis stage.

im0 = im0 + (mx0-mn0)/1000*randn(size(im0));

%% Build the steerable pyramid
[pyr0,pind0] = buildSCFpyr(im0,Nsc,Nor-1);

if ( any(vectify(mod(pind0,2))) )
  error('Algorithm will fail: Some bands have odd dimensions!');
end

%% Subtract mean of lowBand:
nband = size(pind0,1);
pyr0(pyrBandIndices(pind0,nband)) = ...
    real(pyrBand(pyr0,pind0,nband)) - mean2(real(pyrBand(pyr0,pind0,nband)));

rpyr0 = real(pyr0);
apyr0 = abs(pyr0);

%figure(gcf)
%clf
%showIm(im0,'auto',1); title('Original');  drawnow

%% Subtract mean of magnitude:
magMeans0 = zeros(size(pind0,1), 1);
for nband = 1:size(pind0,1)
  indices = pyrBandIndices(pind0,nband);
  magMeans0(nband) = mean2(apyr0(indices));
  apyr0(indices) = apyr0(indices) - magMeans0(nband);
end

%% Compute central autoCorr of lowband

% First make space for all auto-correlations
acr = NaN * ones(Na,Na,Nsc+1);

% get how many bands there are
nband = size(pind0,1);
% and extract the low-resolution residiual channel from the pyramid representation of the image
ch = pyrBand(pyr0,pind0,nband);
% not sure what is going on here - why build a pyramid on this and then extract again? 
% Looks like it is taking the real part of the channel and then splitting again into
% high and low frequency bands with the pyramid and then using the low pass band as the base image
[mpyr,mpind] = buildSFpyr(real(ch),0,0);
im = pyrBand(mpyr,mpind,2);
% get size of channel
[Nly Nlx] = size(ch);
% size of the lowest dimension
Sch = min(Nly,Nlx); %size of low band
% Number of auto-correlation coefficients
le = min(Sch/2-1,la);
cy = Nly/2+1;
cx = Nlx/2+1;
% compute auto-correlation in fourier domain
ac = fftshift(real(ifft2(abs(fft2(im)).^2)))/prod(size(ch));
% keep for later display, full auto-correlation and what it was an auto-correlation of
autoCorrelation{Nsc+1} = ac;
autoCorrelationOf{Nsc+1} = im;
% keep only the required part of the auto-correlation (i.e. the amount specified in the input argument Na)
ac = ac(cy-le:cy+le,cx-le:cx+le);
% put into a big matrix of auto correlation coefficients
acr(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+1) = ac;
% initialize skew and kurtosis
skew0p = zeros(Nsc+1,1);
kurt0p = zeros(Nsc+1,1);
% variance of this image based on the residual is taken from auto-correlation
vari = ac(le+1,le+1);
if vari/var0 > 1e-6,
  % compute skew / kurtosis - the image is already mean subtracted
  % so just use the typical formula cube / or 4th of difference from
  % mean divided by the standard deviation to the power 3 or 4
  skew0p(Nsc+1) = mean2(im.^3)/vari^1.5;
  kurt0p(Nsc+1) = mean2(im.^4)/vari^2;
else
  % unstable esitmate of variance, so use default values for a gaussian distribution
  skew0p(Nsc+1) = 0;
  kurt0p(Nsc+1) = 3;
end

% keep the reconstructed image for display later
reconstructedImage{Nsc+1} = im;

%% Compute  central autoCorr of each Mag band, and the autoCorr of the
%% combined (non-oriented) band.
ace = NaN * ones(Na,Na,Nsc,Nor);
for nsc = Nsc:-1:1,
  for nor = 1:Nor,
    % get which channel of the pyramid we are working on
    nband = (nsc-1)*Nor+nor+1;
    ch = pyrBand(apyr0,pind0,nband);
    % get the size of the channel
    [Nly, Nlx] = size(ch);
    % smallest dimension of channel
    Sch = min(Nlx, Nly);
    % number of coefficients of auto-correlation to keep.
    le = min(Sch/2-1,la);
    cx = Nlx/2+1;  %Assumes Nlx even
    cy = Nly/2+1;
    % Take auto-correlation, using fourier transform - i.e. multiply
    % the fourier transform of the image by itself and then convert back
    % to the image domain using in the inverse fourier transform
    ac = fftshift(real(ifft2(abs(fft2(ch)).^2)))/prod(size(ch));
    % keep for later display, full auto-correlation and what it was an auto-correlation of
    autoCorrelationOrient{Nor*(nsc-1)+nor} = ac;
    autoCorrelationOrientOf{Nor*(nsc-1)+nor} = ch;
    % and placing into a big matrix
    ac = ac(cy-le:cy+le,cx-le:cx+le);
    ace(la-le+1:la+le+1,la-le+1:la+le+1,nsc,nor) = ac;
  end

  %% Combine orientation bands

  % grab all the orientation bands from the pyramid for this spatial scale
  bandNums = [1:Nor] + (nsc-1)*Nor+1;  %ori bands only
  ind1 = pyrBandIndices(pind0, bandNums(1));
  indN = pyrBandIndices(pind0, bandNums(Nor));
  bandInds = [ind1(1):indN(length(indN))];
  %% Make fake pyramid, containing dummy hi, ori, lo
  % that is, we are creating a "fake pyramid" that includes the response across
  % all orientations for this spatial scale, and defaulting everything else
  % to zero so that we can reconstruct with just this scale's information
  fakePind = [pind0(bandNums(1),:);pind0(bandNums(1):bandNums(Nor)+1,:)];
  fakePyr = [zeros(prod(fakePind(1,:)),1);...
	 rpyr0(bandInds); zeros(prod(fakePind(size(fakePind,1),:)),1);];
  % reconstruct the image using only these orientations
  ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
  % upsample to the same resolution as the image
  im = real(expand(im,2))/4;
  % and add this reconstruction to the image that has already been reconstructed from lower spatial scales
  im = im + ch;  
  % Take auto-correlation, using fourier transform - i.e. multiply
  % the fourier transform of the image by itself and then convert back
  % to the image domain using in the inverse fourier transform
  ac = fftshift(real(ifft2(abs(fft2(im)).^2)))/prod(size(ch));
  % keep for later display, full auto-correlation and what it was an auto-correlation of
  autoCorrelation{nsc} = ac;
  autoCorrelationOf{nsc} = im;
  % keep only necessary part of ac (i.e. that specified by the Na argument)
  ac = ac(cy-le:cy+le,cx-le:cx+le);
  acr(la-le+1:la+le+1,la-le+1:la+le+1,nsc) = ac;
  % variance is the auto-correlation at 0 lag
  vari = ac(le+1,le+1);

  if vari/var0 > 1e-6,
    % calculate skew and kurtosis in usual way - note that vari^1.5, vari^2 is std^3 and std^4 as 
    % required by formulas
    skew0p(nsc) = mean2(im.^3)/vari^1.5;
    kurt0p(nsc) = mean2(im.^4)/vari^2;
  else
    % use gaussian defaults when variance estimate is unstable
    skew0p(nsc) = 0;
    kurt0p(nsc) = 3;
  end
  % keep the reconstructed image for display later
  reconstructedImage{nsc} = im;
end

%% Compute the cross-correlation matrices of the coefficient magnitudes
%% pyramid at the different levels and orientations

% initialize correlation matrices, C0 and Cx0 are for correlation of real components
% where C0 is correlation between orientations and Cx0 is between scales. 
% Not sure, but I think these are getting initialized incorrectly - I think
% C0 should have last dimension Nsc and not Nsc+1 since you do the cross-correlations
% inside each scale (and dnon't do anything on the high-pass residual which doesn't have
% orienation bands).
% Likewise, Cx0 should be Nsc-1, since you can't compute a correlation with a lower spatial
% scale once you get to the last scale. The entries into these parts of the matrix seem
% to stay 0 after the code runs.
C0 = zeros(Nor,Nor,Nsc+1);
Cx0 = zeros(Nor,Nor,Nsc);

% same thing but for real components of filter responses
Cr0 = zeros(2*Nor,2*Nor,Nsc+1);
Crx0 = zeros(2*Nor,2*Nor,Nsc);

% iterate over scales
for nsc = 1:Nsc,
  % Get the pyramid band for the first orientation at this spatial scale
  firstBnum = (nsc-1)*Nor+2;
  % size of each orientation band
  cousinSz = prod(pind0(firstBnum,:));
  % calculate the indexes in the pyramid of all the orientation bands at this scale
  ind = pyrBandIndices(pind0,firstBnum);
  cousinInd = ind(1) + [0:Nor*cousinSz-1];

  % if we have not yet reached the top of the pyramid in spatial scale
  if (nsc<Nsc)
    % initialze 
    parents = zeros(cousinSz,Nor);
    rparents = zeros(cousinSz,Nor*2);
    % iterate over orientations
    for nor=1:Nor,
      % get the next lower scale band, at each orientation
      nband = (nsc-1+1)*Nor+nor+1;
      % expand it to the same size as this correlation size
      tmp = expand(pyrBand(pyr0, pind0, nband),2)/4;
      % break into readl and imaginary parts
      rtmp = real(tmp); itmp = imag(tmp);
      %% Double phase:
      tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
      % and put the real and imaginary parts into matrix (vectify is just a function
      % call for turning into a vector and is equaivalent to the matrix(:) notation
      rparents(:,nor) = vectify(real(tmp));
      rparents(:,Nor+nor) = vectify(imag(tmp));
      % get the absolute value
      tmp = abs(tmp);
      % and store that away
      parents(:,nor) = vectify(tmp - mean2(tmp));
    end
  else
    % here we deal with the lowest level of the pyramid in which you cannot go any further
    % just expand out the lowest band to the correct dimensions
    tmp = real(expand(pyrLow(rpyr0,pind0),2))/4;
    % and store in the right place in the matrix
    rparents = [vectify(tmp),...
		vectify(shift(tmp,[0 1])), vectify(shift(tmp,[0 -1])), ...
		vectify(shift(tmp,[1 0])), vectify(shift(tmp,[-1 0]))];
    parents = [];
  end

  % now that everything is in place, compute cross-correlations
  % the terminology here is that cousins are same scale different orientations
  % and parents are across different scales - I believe the parent is
  % the next lower scale, which seems a bit counter-intuitive to me.

  % grab the data for all the orientation bands at this scale
  cousins = reshape(apyr0(cousinInd), [cousinSz Nor]);
  nc = size(cousins,2);   np = size(parents,2);
  % and compute correlation 
  C0(1:nc,1:nc,nsc) = innerProd(cousins)/cousinSz;
  % if there are two scales available, then compute the correlations
  % between scales
  if (np > 0)
    % compute the cross-correlation between scales
    Cx0(1:nc,1:np,nsc) = (cousins'*parents)/cousinSz;
    % at the bottom of the pyramid
    if (nsc==Nsc)
      % just compute the correlation of the parents with itself
      % does this ever get called?
      C0(1:np,1:np,Nsc+1) = innerProd(parents)/(cousinSz/4);
      keyboard
    end
  end
  
  % now do the same thing for the real parts
  cousins = reshape(real(pyr0(cousinInd)), [cousinSz Nor]);
  nrc = size(cousins,2);   nrp = size(rparents,2);  
  Cr0(1:nrc,1:nrc,nsc) = innerProd(cousins)/cousinSz;
  if (nrp > 0)
    Crx0(1:nrc,1:nrp,nsc) = (cousins'*rparents)/cousinSz;
    if (nsc==Nsc)
      Cr0(1:nrp,1:nrp,Nsc+1) = innerProd(rparents)/(cousinSz/4);
    end
  end
end

%% Calculate the mean, range and variance of the LF and HF residuals' energy.

channel = pyr0(pyrBandIndices(pind0,1));
vHPR0 = mean2(channel.^2);

statsLPim = [skew0p kurt0p];

params = struct('pixelStats', statg0, ...
                'pixelLPStats', statsLPim, ...
                'autoCorrReal', acr, ...
                'autoCorrMag', ace, ...
		'magMeans', magMeans0, ...
                'cousinMagCorr', C0, ...
                'parentMagCorr', Cx0, ...
		'cousinRealCorr', Cr0, ...
		'parentRealCorr', Crx0, ...
		'varianceHPR', vHPR0);

% keep the images that we will display in dispPSStats
intermediateSteps.reconstructedImage = reconstructedImage;
intermediateSteps.autoCorrelation = autoCorrelation;
intermediateSteps.autoCorrelationOf = autoCorrelationOf;
intermediateSteps.autoCorrelationOrient = autoCorrelationOrient;
intermediateSteps.autoCorrelationOrientOf = autoCorrelationOrientOf;
