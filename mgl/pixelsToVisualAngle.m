function output = pixelsToVisualAngle(inputs, displaySize)
%%% convert pixels to visual angle
%%% 2022, jiwon yeon
%%% usage
%%%     - pixelsToVisualAngle(inputs, displaySize)
%%% inputs 
%%%     - [x, y] coordiantes or a vector in pixels
%%% displaySize
%%%     - For coordinates. Default displaySize is the screen size

if isvector(inputs)
    if size(inputs,1) == 2 
        inputs = inputs';
        needTranspose = 1;
    elseif all(size(inputs) > 2) || length(size(inputs)) > 2
        error('Inputs must be a scalar or 2D coordinates')
    else
        needTranspose = 0;
    end
else
    needTranspose = 0;
end

if nargin < 2
    displaySize = [mglGetParam('screenWidth'), mglGetParam('screenHeight')];
end
if isscalar(displaySize)
    displaySize = [displaySize, displaySize];
end

xPix2Deg = mglGetParam('xPixelsToDevice');
yPix2Deg = mglGetParam('yPixelsToDevice');

% convert pixels to cm
if isscalar(inputs)
    mean_PixelsInCm = mean([xPix2Deg, yPix2Deg]);
    output = inputs .* mean_PixelsInCm;
else
    % move the centers first
    inputs = [inputs(:,1) - ceil(displaySize(1)/2), ...
        inputs(:,2) - ceil(displaySize(2)/2)];
    % convert it to cm
    output = [inputs(:,1) .* yPix2Deg, ...
        inputs(:,2) .* xPix2Deg];
end

% transpose if needed
if needTranspose
    output = output';
end
    


