%% FSL Combine pe0 and pe1 calibration scans (mux8 acquisitions)
% Based on code from Bob Dougherty (CNI)
% Dan Birman (2015-05)
% dbirman@stanford.edu
%
% Call: [str, unwarp] = fsl_pe0pe1('/path/to/your/directory',doUnwarp=False)
%       returns only a string which can be disp() to see what will be
%       unwarped.
%
% Call: fsl_pe0pe1('/path/to/your/directory',doUnwarp=True,unwarp)
%       Performs the unwarping for the files in 'unwarp' (from a previous
%       doUnwarp=False call).
%
% Failure: str == 'failure'
%
% CODE FROM BOB:
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

function [str, unwarp] = fsl_pe0pe1(folder,unwarp)


[s, r] = system('fslroi');
if s==127
    str = 'failed';
    disp('(fsl_pe0pe1) FSL may not be properly installed. Check your PATH');
    return
end

if ~exist('unwarp')
    doUnwarp = 0;
else
    doUnwarp = 1;
end

if doUnwarp
    disp(sprintf('(fsl_pe0pe1) Unwarping in %s',folder));
    tic
else
    disp(sprintf('(fsl_pe0pe1) Checking %s for unwarp files',folder));
end

if ~exist(folder)
    disp('(fsl_pe0pe1) Folder doesn''t exist... Failed');
    str = 'failed';
    return
end

tfolder = fullfile(folder,'temp');
mkdir(tfolder);
files = dir(folder);

findFiles = 0;
if ~exist('unwarp')
    unwarp.calfiles = {};
    unwarp.EPIfiles = {};
    findFiles = 1;
end

% Currently I wrote this to find one pe1 file, and use that for all of the
% pe0 files.

found_acq_params = 0;

%% Separate files by type
if findFiles
    for i = 1:length(files)
        fi = files(i);
        if strfind(fi.name,'acq_params.txt')
            found_acq_params = 1;
        elseif strfind(fi.name,'uw_')
            % skip files that might have already been unwarped
        elseif ~isempty(strfind(fi.name,'pe1')) || ~isempty(strfind(fi.name,'CAL'))
            unwarp.EPIfiles{end+1} = fi.name;
        elseif ~isempty(strfind(fi.name,'pe0')) || ~isempty(strfind(fi.name,'mux8'))
            if fi.bytes > 100000000 % 1 mega byte
                unwarp.calfiles{end+1} = fi.name;
            else
                disp(sprintf('(fsl_pe0pe1) File size < 100 mB. Likely a cancelled scan: %s',fi.name));
                if strcmp(input('Include? [y/n]: ','s'),'y')
                    unwarp.calfiles{end+1} = fi.name;
                end
            end
        end
    end

    %% Add acq_params.txt if not found

    if ~found_acq_params
        acqFile = fullfile(folder,'acq_params.txt');
        system(sprintf('echo ''0 1 0 1'' > %s',acqFile));
        system(sprintf('echo ''0 -1 0 1'' >> %s',acqFile));
    end

    %% Check if we have multiple pe1 scans

    if length(unwarp.EPIfiles) > 1
        disp('(fsl_pe0pe1) Multiple pe1 scans found, using last scan.\nYou can implement different functionality...');
        scanchoice = 0;
        while ~scanchoice
            in = input('(fsl_pe0pe1) Use first or last scan? [f/l]','s');
            if strcmp(in,'f')
                unwarp.EPIfiles = unwarp.EPIfiles(1); scanchoice=1;
            elseif strcmp(in,'l')
                unwarp.EPIfiles = unwarp.EPIfiles(end); scanchoice=1;
            else
                in = input('(fsl_pe0pe1) Incorrect input. Use first or last scan? [f/l]','s');
            end
        end
    end
else

    acqFile = fullfile(folder,'acq_params.txt');
    if ~isfile(acqFile)
        system(sprintf('echo ''0 1 0 1'' > %s',acqFile));
        system(sprintf('echo ''0 -1 0 1'' >> %s',acqFile));
    end
end

%% Check whether to continue

str = '';
str = strcat(str,'\n','**********************************************************');
for i = 1:length(unwarp.calfiles)
    str = strcat(str,'\n',sprintf('Unwarping %s',unwarp.calfiles{i}));
end
str = strcat(str,'\n',sprintf('Using calibration file %s',unwarp.EPIfiles{1}));
str = strcat(str,'\n','**********************************************************');
str = sprintf(str);

if ~doUnwarp
    return
end

% If we get here we are unwarping!
disp(str);

%% Run fslroi

% pe1
roi1files = {};
roi1files{1} = hlpr_fslroi(unwarp.EPIfiles{1},1,1,1,1,tfolder,folder);
% pe0
roi0files = {};
disppercent(-inf,'Calculating ROIs...');
for i = 1:length(unwarp.calfiles)
    roi0files{i} = hlpr_fslroi(unwarp.calfiles{i},i,0,1,1,tfolder,folder);
    disppercent(i/length(unwarp.calfiles));
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
    finalfiles{i} = hlpr_applytopup(tufiles{i},unwarp.calfiles{i},tfolder,folder);
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

%% Backup + rename
% Move all of the original files (unwarp.calfiles) to folder//unwarp_orig
% Then renamne the uw_ files to the original names
if ~isdir(fullfile(folder,'unwarp_orig'))
    mkdir(fullfile(folder,'unwarp_orig'));
end

for i = 1:length(unwarp.calfiles)
    file = unwarp.calfiles{i};
    fileLoc = fullfile(folder,file);
    backupLoc = fullfile(folder,'unwarp_orig',file);
    % make backup
    system(sprintf('mv %s %s',fileLoc,backupLoc));
    % uw_ file
    uw_fileLoc = fullfile(folder,strcat('uw_',file));
    % rename uw_ file
    system(sprintf('mv %s %s',uw_fileLoc,fileLoc));
end

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