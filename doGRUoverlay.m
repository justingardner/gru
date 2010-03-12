
function doGRUoverlay
%% Convert venogram to an overlay (intensities range from 0-1)
d=cbiReadNifti('veno.img');
vmax=max(d(:));
vmin=min(d(:));
d1=d-vmin;
dif=vmax-vmin;
over=d1/dif;
hdr=cbiReadNiftiHeader('veno.hdr');
over=cbiWriteNifti('venoOVER.img',over,'hdr');

%% edit header to spoof for mrTools

hdr.dim(1,1)=4
hdr.dim(5,1)=1
hdr=cbiWriteNiftiHeader(hdr,'venoOVER.hdr')


%% making overlay directory and moving files
mkdir('Raw/TSeries')
movefile('venoOVER.img','Raw/TSeries/venoOVER.img')
movefile('venoOVER.hdr','Raw/TSeries/venoOVER.hdr')

%% run mrinit
mrInit

%% load overlay in mrLoadRet
disp=cbiReadNifti('Raw/TSeries/venoOVER.img');
mrDispOverlay(disp,1,1);
