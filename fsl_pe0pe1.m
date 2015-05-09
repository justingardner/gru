%% FSL Combine pe0 and pe1 calibration scans (mux8 acquisitions)
% Based on code from Bob Dougherty (CNI)
% Dan Birman (2015-05)
% dbirman@stanford.edu
%
% Call: fsl_pe0pe1('/path/to/your/directory')
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

disp(sprintf('(fsl_pe0pe1) Unwarping in %s',folder));
tic

tfolder = fullfile(folder,'temp');
mkdir(tfolder);
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
    elseif strfind(fi.name,'uw_')
        % skip
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
    system(sprintf('echo ''0 1 0 1'' > %s',acqFile));
    system(sprintf('echo ''0 -1 0 1'' >> %s',acqFile));
end

%% Check if we have multiple pe1 scans

if length(pe1files) > 1
    disp('(fsl_pe0pe1) Multiple pe1 scans found, using first scan.\nYou can implement different functionality...');
    pe1files = pe1files(1);
end

%% Check whether to continue

disp('**********************************************************');
for i = 1:length(pe0files)
    disp(sprintf('Unwarping %s',pe0files{i}));
end
disp(sprintf('Using calibration file %s',pe1files{1}));
disp('**********************************************************');
disp('\n\nIs this correct? [y/n] ');
            
if ~strcmp(input('Include? [y/n]: ','s'),'y')
    disp('(fsl_pe0pe1) Canceling...');
    return
end

%% Run fslroi

% pe1
roi1files = {};
roi1files{1} = hlpr_fslroi(pe1files{1},1,1,1,1,tfolder,folder);
% pe0
roi0files = {};
disppercent(-inf,'Calculating ROIs...');
for i = 1:length(pe0files)
    roi0files{i} = hlpr_fslroi(pe0files{i},i,0,1,1,tfolder,folder);
    disppercent(i/length(pe0files));
end
disppercent(inf);

%% fslmerge
mergefiles = {};
disppercent(-inf,'Merging pe0 and pe1 files...');
for i = 1:length(roi0files)
    mergefiles{i} = hlpr_fslmerge(roi0files{i},roi1files{1},i,tfolder);
    disppercent(i/length(roi0files));
end
disppercent(inf);

%% topup
tufiles = {};
disppercent(-inf,'Calculating topup...');
for i = 1:length(mergefiles)
    tufiles{i} = hlpr_topup(mergefiles{i},i,tfolder,folder);    
    disppercent(i/length(mergefiles));

end
disppercent(inf);

%% applytopup
disppercent(-inf,'Applying topup...');
finalfiles = {};
for i = 1:length(tufiles)
    finalfiles{i} = hlpr_applytopup(tufiles{i},pe0files{i},tfolder,folder);
    disppercent(i/length(tufiles));
end
disppercent(inf);


%% gunzip
disppercent(-inf,'Unzipping and removing .gz files...');
failed = [];
for i = 1:length(finalfiles)
    fi = finalfiles{i};
    [s,r] = system(sprintf('gunzip %s.nii.gz',fi));
    if strfind(r,'No such file')
        disp(sprintf('Unzip failed for file %s',fi));
        failed = [failed i];
    end
    disppercent(i/length(finalfiles));
end
disppercent(inf);

%% rm

if ~isempty(failed)
    disp('Some files failed to unzip, check before we remove anything...');
    keyboard
end

for i = 1:length(finalfiles)
    fi = finalfiles{i};
    system(sprintf('rm %s.nii.gz',fi));
end

%% cleanup
% for i = 1:length(roi1files)
%     system(sprintf('rm %s.nii.gz',fullfile(folder,roi1files{i})));
% end
% for i = 1:length(mergefiles)
%     system(sprintf('rm %s.nii.gz',fullfile(folder,mergefiles{i})));
% end
system(sprintf('rm %s',fullfile(folder,'acq_params.txt')));
system(sprintf('rm -rf %s',tfolder));

%% disp result
T = toc;
disp(sprintf('(fsl_pe0pe1) Unwarping completed successfully for %s',folder));
disp(sprintf('(fsl_pe0pe1) Elapsed time %04.2f s',T));

function outfile = hlpr_applytopup(tu,orig,tfolder,folder)
% applytopup --imain=rs_pe0 --inindex=1 --method=jac --datain=acq_param.txt --topup=rs_topup --out=rs0_unwarped

outfile = fullfile(folder,sprintf('uw_%s',orig(1:end-4)));
tu = fullfile(tfolder,tu);
acqFile = fullfile(folder,'acq_params.txt');
command = sprintf('applytopup --imain=%s --inindex=1 --method=jac --datain=%s --topup=%s --out=%s',fullfile(folder,orig(1:end-4)),acqFile,tu,outfile);
system(command);

function outfile = hlpr_topup(merge,pos,tfolder,folder)

outfile = sprintf('topup_%02.0f',pos);
outfull = fullfile(tfolder,outfile);
merge = fullfile(tfolder,merge);
acqFile = fullfile(folder,'acq_params.txt');
command = sprintf('topup --imain=%s --datain=%s --config=b02b0.cnf --out=%s',merge,acqFile,outfull);
system(command);

function outfile = hlpr_fslmerge(scan0,scan1,pos,folder)

outfile = sprintf('merge_%02.0f',pos);
outfull = fullfile(folder,outfile);
scan0 = fullfile(folder,scan0);
scan1 = fullfile(folder,scan1);
command = sprintf('fslmerge -t %s %s %s',outfull,scan0,scan1);
system(command);

function outfile = hlpr_fslroi(scan,pos,type,n1,n2,tfolder,folder)

outfile = sprintf('pe%i_%02.0f',type,pos);
outfull = fullfile(tfolder,outfile);
scan = fullfile(folder,scan);
command = sprintf('fslroi %s %s %i %i',scan,outfull,n1,n2);
system(command);