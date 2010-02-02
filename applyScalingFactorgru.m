function applyScalingFactor(baseImg,sFacs)
%rename old file with extension preGradCoil
moveoldfilenametothis = setext(sprintf('%s_preGradCoil',stripext(baseImg)),'hdr');

%check if this function has been run in the past
if isfile(moveoldfilenametothis);
    disp(sprintf('Pre gradient coil replacement header file found in this directory. You may have already scaled this image. Image not scaled.'));
    return
end

sMat = diag(sFacs);
sMat = [[sMat; 0 0 0;] [0 0 0 1]'];

basehdr = setext(baseImg,'hdr');
[data hdr] = cbiReadNifti(basehdr);

newsform44 = sMat*hdr.sform44;
newqform44 = sMat*hdr.qform44;
newhdr = cbiSetNiftiSform(hdr,newsform44);
newhdr = cbiSetNiftiQform(hdr,newqform44);

%move the old header to have this extention
command = sprintf('!mv %s %s',baseImg,moveoldfilenametothis);
eval(command);

if isfile(moveoldfilenametothis)
    newfilename = setext(sprintf('%s',stripext(baseImg)),'hdr');
    cbiWriteNifti(newfilename,data,newhdr);
end

end