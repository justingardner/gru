% dicom2xform.m
%
%        $Id:$ 
%      usage: xform = dicom2xform(dicomHeader) 
%         by: justin gardner
%       date: 05/21/15
%    purpose: Gets a nifti-style xform from dicomHeader, based on code from Bob Dougherty
%             dicomHeader can be a header returned by dicominfo or a filename of a dicom
%
function xform = dicom2xform(dicomHeader)

xform = [];

% try to load from file if string
if isstr(dicomHeader)
  if isfile(dicomHeader)
    dicomHeader = dicominfo(dicomHeader);
  else
    disp(sprintf('(dicom2xform) No file found: %s',dicomHeader));
    return
  end
end

% fields to check and how long they should be
checkFields = {'PixelSpacing','SpacingBetweenSlices','ImagePositionPatient','ImageOrientationPatient'};
fieldLengths = [2 1 3 6];
for iField = 1:length(checkFields)
  if ~isfield(dicomHeader,checkFields{iField})
    disp(sprintf('(dicom2xform) Field %s not found',checkFields{iField}));
    return
  elseif length(dicomHeader.(checkFields{iField})) ~= fieldLengths(iField)
    disp(sprintf('(dicom2xform) Field %s expected to have length %i but had %i',fieldLengths(iField),length(dicomHeader.(checkFields{iField}))));
    return
  end
end

% ok get spacing
xSpacing = dicomHeader.PixelSpacing(1);
ySpacing = dicomHeader.PixelSpacing(2);

% get spacing between slices
zSpacing = dicomHeader.SpacingBetweenSlices;

% get voxel dimensions (these are in mm)
voxDim = [xSpacing ySpacing zSpacing];

% get origin
origin = dicomHeader.ImagePositionPatient(:) .* [-1 -1 1]';

% get orientation variables
rowCos = dicomHeader.ImageOrientationPatient(1:3);
colCos = dicomHeader.ImageOrientationPatient(4:6);
sliceNorm = cross(rowCos,colCos);

% compute rotationMatrix
rotMatrix = [-rowCos(1) -colCos(1) -sliceNorm(1);...
	     -rowCos(2) -colCos(2) -sliceNorm(2);...
	      rowCos(3)  colCos(3)  sliceNorm(3)];

% now create the xform
xform = zeros(4,4);
xform(1:3,1:3) = rotMatrix*diag(voxDim);
xform(1:3,4) = origin;
xform(4,4) = 1;

