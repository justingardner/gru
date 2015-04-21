%% FSL Combine pe0 and pe1 calibration scans (mux8 acquisitions)
% Base on code from Bob Dougherty (CNI)
% Call: fsl_pe0pe1(directory)
%
% #!/bin/bash
% 
% # topup for rs data mux8 hcp resting state data:
% 
% # To compute the echo train length, run:
% # fslhd rs_pe0.nii.gz | grep desc
% # and compute acq[0]*ec/1000
% echo '0 1 0 0.05720' > acq_params.txt
% echo '0 -1 0 0.05720' >> acq_params.txt
% fslroi rs_pe0.nii.gz bu 1 1
% fslroi rs_pe1.nii.gz bd 1 1
% fslmerge -t bud bu bd
% topup --imain=bud --datain=acq_param.txt --config=b02b0.cnf --out=rs_topup
% applytopup --imain=rs_pe0 --inindex=1 --method=jac --datain=acq_param.txt --topup=rs_topup --out=rs0_unwarped
% applytopup --imain=rs_pe1 --inindex=2 --method=jac --datain=acq_param.txt --topup=rs_topup --out=rs1_unwarped
%
% FSL Requests that we include the following text in any manuscripts that
% use this function:
% Brief summary text: "Data was collected with reversed phase-encode blips, resulting in pairs of images with distortions going in opposite directions. From these pairs the susceptibility-induced off-resonance field was estimated using a method similar to that described in [Andersson 2003] as implemented in FSL [Smith 2004] and the two images were combined into a single corrected one."
% 
% [Andersson 2003] J.L.R. Andersson, S. Skare, J. Ashburner How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 20(2):870-888, 2003.
% 
% [Smith 2004] S.M. Smith, M. Jenkinson, M.W. Woolrich, C.F. Beckmann, T.E.J. Behrens, H. Johansen-Berg, P.R. Bannister, M. De Luca, I. Drobnjak, D.E. Flitney, R. Niazy, J. Saunders, J. Vickers, Y. Zhang, N. De Stefano, J.M. Brady, and P.M. Matthews. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 23(S1):208-219, 2004. 

function fsl_pe0pe1(folder)

% disp(sprintf('(fsl_pe0pe1) Warning: This function takes a LONG time to run locally.\nRunning this on the LXC server would be much better and could be run in parallel.'));

files = dir(folder);

pe0files = {};
pe1files = {};

% Currently I wrote this to find one pe1 file, and use that for all of the
% pe0 files.

found_acq_params = 0;

%% Separate files by type
for i = 1:length(files)
    fi = files(i);
    if strfind(fi.name,'acq_params.txt')
        found_acq_params = 1;
    elseif strfind(fi.name,'pe0')
        if fi.bytes > 100000000 % 1 mega byte
            pe0files{end+1} = fi.name;
        else
            disp(sprintf('(fsl_pe0pe1) File size < 100 mB. Likely a cancelled scan: %s',fi.name));
            if strcmp(input('Include? [y/n]: ','s'),'y')
                pe0files{end+1} = fi.name;
            end
        end
    elseif strfind(fi.name,'pe1')
        pe1files{end+1} = fi.name;
    end
end

%% Add acq_params.txt if not found

if ~found_acq_params
    acqFile = fullfile(folder,'acq_params.txt');
    system(sprintf('echo ''0 1 0 0.05720'' > %s',acqFile));
    system(sprintf('echo ''0 -1 0 0.05720'' >> %s',acqFile));
end

%% Check if we have multiple pe1 scans

if length(pe1files) > 1
    disp('(fsl_pe0pe1) Multiple pe1 scans found, using first scan.\nYou can implement different functionality...');
    pe1files = pe1files(1);
end
%% Run fslroi

% pe1
roi1files = {};
roi1files{1} = hlpr_fslroi(pe1files{1},1,1,1,1,folder);
% pe0
roi0files = {};
for i = 1:length(pe0files)
    roi0files{i} = hlpr_fslroi(pe0files{i},i,0,1,1,folder);
end

%% fslmerge
mergefiles = {};
for i = 1:length(roi0files)
    mergefiles{i} = hlpr_fslmerge(roi0files{i},roi1files{1},i,folder);
end

%% topup
tufiles = {};
for i = 1:length(mergefiles)
    tufiles{i} = hlpr_topup(mergefiles{i},i,folder);
end

%% applytopup
finalfiles = {};
for i = 1:length(tufiles)
    finalfiles{i} = hlpr_applytopup(tufiles{i},pe0files{i},folder);
end

%% cleanup
for i = 1:length(roi1files)
    system(sprintf('rm %s',fullfile(folder,roi1files{i})));
end
for i = 1:length(mergefiles)
    system(sprintf('rm %s',fullfile(folder,mergefiles{i})));
end
system(sprintf('rm %s',fullfile(folder,'acq_params.txt')));

%% disp result
disp(sprintf('(fsl_pe0pe1) Completed successfully for %s',folder));

function outfile = hlpr_applytopup(tu,orig,folder)
% applytopup --imain=rs_pe0 --inindex=1 --method=jac --datain=acq_param.txt --topup=rs_topup --out=rs0_unwarped

outfile = fullfile(folder,sprintf('uw_%s',orig));
tu = fullfile(folder,tu);
command = sprintf('applytopup --imain=%s --inindex=1 --method=jac --datain=acq_params.txt --topup=%s --out=%s',orig,tu,outfile);
system(command);

function outfile = hlpr_topup(merge,pos,folder)

outfile = sprintf('topup_%02.0f',pos);
outfull = fullfile(folder,outfile);
merge = fullfile(folder,merge);
command = sprintf('topup --imain=%s --datain=acq_params.txt --config=b02b0.cnf --out=%s',merge,outfull);
system(command);

function outfile = hlpr_fslmerge(scan0,scan1,pos,folder)

outfile = sprintf('merge_%02.0f',pos);
outfull = fullfile(folder,outfile);
scan0 = fullfile(folder,scan0);
scan1 = fullfile(folder,scan1);
command = sprintf('fslmerge -t %s %s %s',outfull,scan0,scan1);
system(command);

function outfile = hlpr_fslroi(scan,pos,type,n1,n2,folder)

outfile = sprintf('pe%i_%02.0f',type,pos);
outfull = fullfile(folder,outfile);
scan = fullfile(folder,scan);
command = sprintf('fslroi %s %s %i %i',scan,outfull,n1,n2);
system(command);