function output = pixelsToVisualAngle(inputs, displaySize)
%%% convert pixels to visual angle
%%% 2022, jiwon yeon
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

x_PixelsInCm = mglGetParam('xPixelsToDevice');
y_PixelsInCm = mglGetParam('yPixelsToDevice');
screen_distance = mglGetParam('devicePhysicalDistance');

% convert pixels to cm
if isscalar(inputs)
    mean_PixelsInCm = mean([x_PixelsInCm, y_PixelsInCm]);
    inputs_cm = inputs .* mean_PixelsInCm;
else
    % move the centers first
    inputs = [inputs(:,1) - ceil(displaySize(1)/2), ...
        inputs(:,2) - ceil(displaySize(2)/2)];
    % convert it to cm
    inputs_cm = [inputs(:,1) .* x_PixelsInCm, ...
        inputs(:,2) .* y_PixelsInCm];
end

% convert cm to visual angle
output = atan(inputs_cm ./ (2 * screen_distance)) .* 2 .* (180/pi);

% transpose if needed
if needTranspose
    output = output';
end
    


