function [sFacmm sFacvox] = findScalingFactor(filename,useQform)

if ieNotDefined('useQform'), useQform = 0; end

filenamehdr = setext(filename,'hdr');
hdr = cbiReadNiftiHeader(filenamehdr);

sinv = inv(hdr.sform44);
qinv = inv(hdr.qform44);
pixdim = hdr.pixdim(2:4);

%find out what a shift of 1 in the magnet means for each axis of old image
Amag = [0 0 0 1]'; %position 0,0,0 in magnet
%BmagmmBig = [1 0 0; 0 1 0; 0 0 1]; %positions shifted 1mm in each axis
%BmagBig = diag(BmagmmBig*(1./pixdim)); %the number of voxels equal to a shift of 1mm
%BmagBig = [BmagBig [1 1 1]']; %dimension adjustment for matrix multiplication
BmagBig = [1 0 0 1; 0 1 0 1; 0 0 1 1]; %positions shifted '1' in each axis

for i = 1:3 %for each axis (each row in BmagBig)
    
    Bmag = BmagBig(i,:)'; %column vector 1001,0101,0011 in mm
    if useQform == 1, xforminv = qinv; else  xforminv = sinv; end
    Aimg = xforminv*Amag; %what no shift means for this image 
    Bimg = xforminv*Bmag; %what a shift of 1 means for this image
    Dimg = Aimg-Bimg; %what the raw distance of a shift of 1mm outside magnet means for the image
    Dimg = Dimg(1:3); %shed the 4th row - it's a zero
    Distvox(i) = sqrt(sum(Dimg.^2)); %this is the image distance in mm for a shift of 1mm outside the magnet
    
    %multiply the raw distance by the pixdim to get this in mm...
    Distmm(i) = sqrt(sum((Dimg.*pixdim(1:3)).^2)); %this is the image distance in vox for a shift of 1 outside the magnet
end
%we now have Distvox = [xvox yvox zvox] and Distmm = [xmm ymm zmm]
sFacmm = Distmm;
sFacvox = Distvox;