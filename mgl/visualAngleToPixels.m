function output = visualAngleToPixels(inputs, displaySize)
%%% convert visual angle to pixels
%%% 2022, jiwon yeon
%%% usage
%%%     - visualangleToPiexels(inputs, displaySize)
%%% inputs 
%%%     - [x, y] coordiantes or a vector in visual angle
%%% displaySize
%%%     - For coordinates, the default display size is the screen size

if nargin < 2
    displaySize = [mglGetParam('screenWidth'), mglGetParam('screenHeight')];
end
if isscalar(displaySize)
    displaySize = [displaySize, displaySize];
end

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

% get number of pixels per 1 degree visual angle on the display
xDeg2Pix = mglGetParam('xDeviceToPixels');
yDeg2Pix = mglGetParam('yDeviceToPixels');

% change visual angle values to pixels
if isscalar(inputs)
    output = round(inputs .* max([xDeg2Pix,yDeg2Pix]));        
else    
    output = [inputs(:,1) .* xDeg2Pix, ...
        inputs(:,2) .* yDeg2Pix];    
    % centering
    output = round([output(:,1) + floor(displaySize(1)/2), ...
        output(:,2) + floor(displaySize(2)/2)]);     
end

if needTranspose
    output = output';
end


