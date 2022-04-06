function output = visualAngleToPixels(inputs, displaySize)
%%% convert visual angle to pixels
%%% 2022, jiwon yeon
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

% get number of pixels per cm on the display
x_nPixPerCm = mglGetParam('xDeviceToPixels');
y_nPixPerCm = mglGetParam('yDeviceToPixels');
screen_distance = mglGetParam('devicePhysicalDistance');

% convert visual angle to cm
inputs_cm = 2 .* screen_distance .* tan(inputs ./2 .* (pi/180));

% change the cm values to pixels
if isscalar(inputs)
    nPixPerCm = mean([x_nPixPerCm, y_nPixPerCm]);
    output = round(inputs .* nPixPerCm);        
else    
    output = round([inputs_cm(:,1) .* x_nPixPerCm, ...
        inputs_cm(:,2) .* y_nPixPerCm]);    
    % centering
    output = [output(:,1) + ceil(displaySize(1)/2), ...
        output(:,2) + ceil(displaySize(2)/2)];     
end

if needTranspose
    output = output';
end


