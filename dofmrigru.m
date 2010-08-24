% dofmrigru.m
%
%        $Id:$ 
%      usage: dofmrigru
%         by: justin gardner
%       date: 07/06/09
%    purpose: do initial processing stream from fid files collected on the Varian magnet
%             to MLR nifti file format. You will run this program twice to
%             process this data. The first time files will get moved into
%             your local computer. You then have to manually generate a
%             mask and run peak. Then you run it a second time to do the
%             sense reconstruction and run physiofix, etc.
%
%             To use this, first make a directory where you want
%             all the data saved on your local machine. Use the same directory name
%             as the one on /usr1 or /usr4 which contains the raw data. For example:
%
%             cd ~/data
%             mkdir S00120091216
%             cd S00120091216
%
%             Now, make sure that your raw data is setup as follows:
%
%             /usr1/justin/data/S00120091216 should contain all of your data. Including
%             all of the .fid directories)
%             all of the car/ext files
%             all of the stimfiles (if you have them)
%             any pdf files (documents if you have them)
% 
%             Now you can run dofmrigru as follows:
%
%             dofmrigru('dataDir=/usr1/justin/data')
%
%             Other options:
%                'epirri=eprri5': set which epirri processing function to use
%                'postproc=pp': set which postproc program to use
%                'sense=sense_mac_intel': set which sense reconstruction to use
%                'fidDir=/usr1/justin/data/s00620101001/Pre': Set this if you want to load the first pass
%                    fid files form a specific directory.
%                'carextDir=/usr1/justin/data/s00620101001/carext': Set this if you want to load the first pass
%                    car/ext files form a specific directory.
%                'stimfileDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    stimfiles form a specific directory.
%                'pdfDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    pdf files form a specific directory.
%                'numMotionComp=1': Set to 0 if you don't want to run MLR motion comp. Set to > 1 if you want
%                    to set multiple motionComp parameters (e.g. for motionComping two sets of scans taken at
%                    different resolutions)
%
%             First pass will sort through the specified datadir and copy
%             all of these files in the correct directories on your local
%             computer.
%
function retval = dofmrigru(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14])
  help dofmrigru
  return
end

% Default arguments
global dataDir;
fidDir = [];
carextDir = [];
pdfDir = [];
stimfileDir = [];
numMotionComp = [];
movepro=[];
global epirri;
global postproc;
global sense;
getArgs(varargin,{'dataDir=/usr1/justin/data','fidDir=[]','carextDir=[]','pdfDir=[]','stimfileDir=[]','epirri=epirri5','postproc=pp','sense=/usr1/mauro/SenseProj/command_line/current/executables/sense_mac_intel','numMotionComp=1','movepro=0'});

% check to make sure we have the computer setup correctly to run epirri, postproc and sense
if checkfMRISupportUnitCommands == 0,return,end

% make sure that MLR is not running
mrQuit;

% see if this is the first preprocessing or the second one
if ~isdir('Pre')
  disp(sprintf('(dofmrigru) Running Initial dofmrigru process'));
  dofmrigru1(fidDir,carextDir,stimfileDir,pdfDir,numMotionComp);
else
  disp(sprintf('(dofmrigru) Running Second dofmrigru process'));
  dofmrigru2(movepro);
end

%%%%%%%%%%%%%%%%%%%%
%%   dofmrigru2   %%
%%%%%%%%%%%%%%%%%%%%
function dofmrigru2(movepro)

% get fid file list
fiddir = 'Pre';
fidList = getFileList(fiddir,'fid');
fidList = getFidInfo(fidList);
fidList = sortFidList(fidList);

% get list of epi scans
epiNums = getEpiScanNums(fidList);

% check for sense processing
senseProcessing = isSenseProcessing(fidList,epiNums);

if senseProcessing
  % look for mask in nifti form
  %the mask is made manually before dofmrigru1 from noise/ref in maskdir.
  %copy the correct maskfile into Pre MANUALLY
  maskList = getFileList(fiddir,'hdr','img','mask'); 
  if isempty(maskList)
    disp(sprintf('(dofrmigru2) Could not find any mask files'));
    return
  end
  noiseList = getFileList(fiddir,'edt','epr','noise');
  if isempty(noiseList)
    disp(sprintf('(dofrmigru2) Could not find any noise files'));
    return
  end
  refList = getFileList(fiddir,'edt','epr','ref');
  if isempty(refList)
    disp(sprintf('(dofrmigru2) Could not find any ref files'));
    return
  end
else
  maskList = [];noiseList = [];refList = [];
end

% get anatomy scan nums
anatNums = getAnatScanNums(fidList);

% check to make sure that all epi runs are processed
[epiNumsWithPeaks epiNumsWithCarExt] = checkForPeaks(fidList,epiNums);

% get some info for sense processing
if senseProcessing
  senseScanNums = getSenseScanNums(fidList,epiNums);
  % get which mask goes for each sense file
  [maskNums noiseNums refNums] = chooseMaskForSenseFiles(fidList,maskList,noiseList,refList,epiNums);
else
  % no sense processing defaults
  senseScanNums = [];maskNums = [];noiseNums = [];refNums = [];
end

% check for something to do
if isempty(epiNumsWithCarExt)
  disp(sprintf('(dofmrigru) No epi scans to process'));
  return
end

% check to see if we should quit given that some of the peak files are missing
if length(epiNumsWithPeaks) ~= length(epiNums)
  if ~askuser('You are missing some peak files. Do you still want to continue'),return,end
end

dispList(fidList,epiNumsWithCarExt,sprintf('Epi scans to process'));
dispList(fidList,anatNums,sprintf('Anatomy files'));

% display the commands for processing
processFiles(1,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro);

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
processFiles(0,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro);

%%%%%%%%%%%%%%%%%%%%%%
%%   processFiles   %%
%%%%%%%%%%%%%%%%%%%%%%
function processFiles(justDisplay,fidList,maskList,noiseList,refList,epiNumsWithPeaks,epiNumsWithCarExt,senseScanNums,anatNums,maskNums,noiseNums,refNums,senseProcessing,movepro)

global postproc;
global sense;

command = sprintf('cd Pre');
if justDisplay,disp(command),else,eval(command),end

% open the logfile
if ~justDisplay
  openLogfile('dofmrigru2.log');
end

% convert anatomy files into nifti
for i = 1:length(anatNums)
  disp(sprintf('=============================================='));
  disp(sprintf('Copy processed anatomy into Anatomy directory'));
  disp(sprintf('=============================================='));
  command = sprintf('copyfile %s ../Anatomy',setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'hdr'));
  if justDisplay,disp(command),else,eval(command),end
  command = sprintf('copyfile %s ../Anatomy',setext(fixBadChars(stripext(fidList{anatNums(i)}.filename),{'.','_'}),'img'));
  if justDisplay,disp(command),else,eval(command),end
end

if senseProcessing
  for i = 1:length(maskList)
    disp(sprintf('=============================================='));
    disp(sprintf('Convert masks from nifti to sdt/spr'));
    disp(sprintf('=============================================='));
    % convert mask from nifti to to sdt/spr form
    maskname = maskList{i}.filename; % still in nifti form
    command = sprintf('[d h] = cbiReadNifti(''%s'');',maskname);
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('writesdt(''%s'',d);',setext(maskname,'sdt'));
    if justDisplay,disp(command),else,eval(command),end
  end
end

% and delete headers
disp(sprintf('=============================================='));
disp(sprintf('Remove temporary nifti files'));
disp(sprintf('=============================================='));
command = 'mysystem(''rm -f ../Raw/TSeries/*.hdr'');';
if justDisplay,disp(command),else,eval(command);end
command = 'mysystem(''rm -f ../Raw/TSeries/*.img'');';
if justDisplay,disp(command),else,eval(command);end

disp(sprintf('=============================================='));
disp(sprintf('Physiofix, sense process and convert to nifti epi files'));
disp(sprintf('=============================================='));
% physiofix the files

for i = 1:length(epiNumsWithCarExt)
  % get some names for things
  fidname = fidList{epiNumsWithCarExt(i)}.filename;
  edtname = setext(fidList{epiNumsWithCarExt(i)}.filename,'edt');
  sdtname = setext(fidList{epiNumsWithCarExt(i)}.filename,'sdt');
  hdrname = setext(fidList{epiNumsWithCarExt(i)}.filename,'hdr');
  imgname = setext(fidList{epiNumsWithCarExt(i)}.filename,'img');

  
  % see if it is sense
  if any(i==senseScanNums)
    % noise, ref and mask numbers
    noiseNum = noiseNums(find(i==senseScanNums));
    refNum = refNums(find(i==senseScanNums));
    maskNum = maskNums(find(i==senseScanNums));
    %covarNum = covarNums(find(i==senseScanNums));
    maskname = setext(maskList{maskNum}.filename,'sdt');
    %covarname = covarList{covarNum}.filename;
    noisename = noiseList{noiseNum}.filename;
    refname = refList{refNum}.filename;
    %check to see if this epi has peaks. if so run postproc with physiofix
    if any(epiNumsWithCarExt(i) == epiNumsWithPeaks)
        % convert to epis to edt file, doing physiofix and dc correction
        command = sprintf('mysystem(''%s -outtype 1 -dc -physiofix %s %s'');',postproc,fidname,edtname);
    else
        % convert to epis to edt file, don't do physiofix but dc correction
        command = sprintf('mysystem(''%s -outtype 1 -dc %s %s'');',postproc,fidname,edtname);
    end
    if justDisplay,disp(command),else,eval(command),end
    % run the sense processing for this file
    command = sprintf('mysystem(''%s -data %s -full %s -noise %s -mask %s -remove -recon %s'');',sense,stripext(fidname),stripext(refname),stripext(noisename),stripext(maskname),stripext(fidname));
    % command = sprintf('mysystem(''%s -? -? %s %s
    % %s'')',sense,maskname,covarname,edtname); this is for super old sense!
    if justDisplay,disp(command),else,eval(command),end
    % then convert the sdt file into a nifti
    command = sprintf('[hdr] = fid2niftihdr(''%s'');',fidname);
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('data = readsdt(''%s.sdt'');',stripext(sdtname));
    if justDisplay,disp(command),else,eval(command),end
    command = sprintf('cbiWriteNifti(''%s.hdr'',data.data,hdr);',stripext(fullfile('..','Raw','TSeries',hdrname)));
    if justDisplay,disp(command),else,eval(command),end
  else
    %check to see if this epi has peaks. if so run postproc with physiofix
    if any(epiNumsWithCarExt(i) == epiNumsWithPeaks)
      % run with physiofix correction
      ppoptions = '-dc -physiofix';
    else
      % run without phsyiofix corrections
      ppoptions = '-dc';
    end
    % if we don't have to movepro, then just use postproc to read
    % the fid directly
    if movepro==0
      % convert epis to sdt file, doing physiofix and dc correction
      command = sprintf('mysystem(''%s -outtype 3 %s %s %s'');',postproc,ppoptions,fidname,imgname);
    else
      % if we have to movepro, then convert the file to sdt using
      % fid2nifti.
      command = sprintf('[d h] = fid2nifti(''%s'',''movepro=%f'');writesdt(''%s'',d);mysystem(''%s -intype 2 -outtype 3 %s %s %s'');',fidname,movepro,setext(fidname,'sdt'),postproc,ppoptions,setext(fidname,'sdt'),imgname);
    end
    if justDisplay,disp(command),else,eval(command),end
    % then convert the output of postproc to a valid nifti file, by
    % pasting on the header from fid2niftihdr
    command = sprintf('[hdr info] = fid2niftihdr(''%s'',1,''movepro=%f'');',fidname,movepro);
    if justDisplay,disp(command),else,eval(command),end
    % add back reference volume to header
    command = sprintf('if info.nRefVolumes,hdr.dim(5) = hdr.dim(5) + info.nRefVolumes;end');
    if justDisplay,disp(command),else,eval(command),end
    % save header
    command = sprintf('cbiWriteNiftiHeader(hdr,''%s'');',hdrname);
    if justDisplay,disp(command),else,eval(command),end
    % now load completed nifti file for saving into Raw/TSeries
    command = sprintf('[d h] = cbiReadNifti(''%s'');',hdrname);
    if justDisplay,disp(command),else,eval(command),end
    % now, remove reference volumes from nifti file if necessary, 
    command = sprintf('if info.nRefVolumes,d = d(:,:,:,info.nRefVolumes+1:end);h.dim(5) = h.dim(5)-info.nRefVolumes;end',hdrname,hdrname);
    if justDisplay,disp(command),else,eval(command),end
    % save nifti to Raw/TSeries
    command = sprintf('cbiWriteNifti(''%s'',d,h);',fullfile('..','Raw','TSeries',hdrname));
    if justDisplay,disp(command),else,eval(command),end
  end
end

command = sprintf('cd ..');
if justDisplay,disp(command),else,eval(command),end

disp(sprintf('=============================================='));
disp(sprintf('Run motion comp'));
disp(sprintf('=============================================='));
if ~justDisplay
  motionCompFilenameNum = 1;
  motionCompFilename = sprintf('motionCompParams1.mat');
  while isfile(motionCompFilename)
    % load the params
    eval(sprintf('load %s',motionCompFilename));
    % run the motion comp
    v = newView;
    v = motionComp(v,params);
    deleteView(v);
    % get the next motionCompParams to run
    motionCompFilenameNum = motionCompFilenameNum+1;
    motionCompFilename = sprintf('motionCompParams%i.mat',motionCompFilenameNum);
  end
end

disp(sprintf('=============================================='));
disp(sprintf('DONE'));
disp(sprintf('=============================================='));

if ~justDisplay
  closeLogfile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   chooseMaskForSenseFiles   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maskNums noiseNums refNums] = chooseMaskForSenseFiles(fidList,maskList,noiseList,refList,epiNums)
% set up variables
maskNums = []; noiseNums = []; refNums = [];
if isempty(epiNums),return,end

% get the names of the mask files
for i = 1:length(maskList)
  maskNames{i} = maskList{i}.filename;
end

% get the names of the noise files
for i = 1:length(noiseList)
  noiseNames{i} = noiseList{i}.filename;
end

% get the names of the covar files
for i = 1:length(refList)
  refNames{i} = refList{i}.filename;
end

% get the names of the covar files
%for i = 1:length(covarList)
%  covarNames{i} = covarList{i}.filename;
%end

% get info from scans
for i = 1:length(epiNums)
  scanNames{i} = stripext(fidList{epiNums(i)}.filename);
  maskNamesMatch{i} = maskNames;
  noiseNamesMatch{i} = noiseNames;
  refNamesMatch{i} = refNames;
  % covarNamesMatch{i} = covarNames;
  startTime{i} = fidList{epiNums(i)}.startTimeStr;
  endTime{i} = fidList{epiNums(i)}.endTimeStr;
end

% make an mrParamsDialog where you can select
paramsInfo{1} = {'scanNum',1,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',length(epiNums))};
paramsInfo{2} = {'scanName',scanNames,'group=scanNum','type=String','editable=0'};
paramsInfo{3} = {'startTime',startTime,'group=scanNum','type=String','editable=0'};
paramsInfo{4} = {'endTime',endTime,'group=scanNum','type=String','editable=0'};
paramsInfo{5} = {'maskFile',maskNamesMatch,'group=scanNum'};
% paramsInfo{6} = {'covarFile',covarNamesMatch,'group=scanNum'};
paramsInfo{6} = {'noiseFile',noiseNamesMatch,'group=scanNum'};
paramsInfo{7} = {'refFile',refNamesMatch,'group=scanNum'};

% put out the dialog
if (length(maskNames) == 1)&&(length(refNames) == 1)&&(length(noiseNames) == 1)% && (length(covarNames) == 1)
  % if there is nothing to choose (i.e. there are not multiple mask and covar, then just select the default)
  params = mrParamsDefault(paramsInfo);
else
  % otherwise have the user choose
  params = mrParamsDialog(paramsInfo,'Select which mask file we want to use');
end
if isempty(params),return;end

% get whcih mask and covar to use
for i = 1:length(epiNums)
  maskNums(i) = find(strcmp(params.maskFile{i},maskNames));
  noiseNums(i) = find(strcmp(params.noiseFile{i},noiseNames));
  refNums(i) = find(strcmp(params.refFile{i},refNames));
  % covarNums(i) = find(strcmp(params.covarFile{i},covarNames));
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   checkForPeaks   %%
%%%%%%%%%%%%%%%%%%%%%%%
function [epiNumsWithPeaks epiNumsWithCarExt] = checkForPeaks(fidList,epiNums)

passedCheckCarExt = ones(1,length(epiNums));
passedCheckPeaks = ones(1,length(epiNums));

for i = 1:length(epiNums)
  % shortcut
  fid = fidList{epiNums(i)};
  % check for car/ext
  if ~all(fid.hasCarExt)
      passedCheckCarExt(i) = 0;
    % display what is missing
    filenames = {'car','ext'};
    for j = 1:length(filenames)
      if ~fid.hasCarExt(j)
	disp(sprintf('%s: Missing %s file',fid.filename,filenames{j}));
      end
    end
  end
  % check for peak files
  if ~all(fid.hasPeakFiles)
    passedCheckPeaks(i) = 0;
    % display what is missing
    filenames = {'acq.peak.bit','cardio.peak.bit','respir.peak.bit'};
    for j = 1:length(filenames)
      if ~fid.hasPeakFiles(j)
	disp(sprintf('%s: Missing %s file (need to run peak)',fid.filename,filenames{j}));
      end
    end
  end
end

scansThatPassedPeaks = find(passedCheckPeaks);
epiNumsWithPeaks = [];
for i = 1:length(scansThatPassedPeaks)
    epiNumsWithPeaks(end+1) = epiNums(scansThatPassedPeaks(i));
end

scansThatPassedCarExt = find(passedCheckCarExt);
epiNumsWithCarExt = [];
for i = 1:length(scansThatPassedCarExt)
    epiNumsWithCarExt(end+1) = epiNums(scansThatPassedCarExt(i));
end

%%%%%%%%%%%%%%%%%%%%
%%   dofmrigru1   %%
%%%%%%%%%%%%%%%%%%%%
function dofmrigru1(fidDir,carextDir,stimfileDir,pdfDir,numMotionComp)

global dataDir;

% get the current directory
expdir = getLastDir(pwd);

% location of all directories. Should be a directory
% in there that has the same name as the current directory
% i.e. s00120090706 that contains the fid files and the
% car/ext files
if isempty(fidDir),fidDir = fullfile(dataDir,expdir);end
if isempty(carextDir),carextDir = fullfile(dataDir,expdir);end
if isempty(stimfileDir), stimfileDir = fullfile(dataDir,expdir);end
if isempty(pdfDir),pdfDir = fullfile(dataDir,expdir);end

% check the fiddir
if ~isdir(fidDir)
  disp(sprintf('(dofmrigru1) Could not find fid directory %s',fidDir));
  return
end

% now get info about fids
[fidList fidListArray] = getFileList(fidDir,'fid');
fidList = getFidInfo(fidList);
fidList = sortFidList(fidList);

% get list of epi scans
epiNums = getEpiScanNums(fidList);
if isempty(epiNums)
  disp(sprintf('(dofmrigru1) Could not find any epi files in %s',fidDir));
  return
end

% check to see if any of the epi scans have acceleration greater than 1
senseProcessing = isSenseProcessing(fidList,epiNums);

% get anatomy scan nums
anatNums = getAnatScanNums(fidList);
if isempty(anatNums)
    disp(sprintf('(dofmrigru1) Could not find any non-raw anatomies in %s',fidDir));
    return
end

% get sense noise/ref scans
[senseNoiseNums senseRefNums] = getSenseNums(fidList);
if isempty(senseRefNums) && senseProcessing
    disp(sprintf('(dofmrigru1) Could not find any sense ref files in %s',fidDir));
    return
end
if isempty(senseNoiseNums) && senseProcessing
    disp(sprintf('(dofmrigru1) Could not find any sense noise files in %s',fidDir));
    return
end


% find car/ext files
carList = getFileList(carextDir,'car','ext');
carList = getCarInfo(carList);
if isempty(carList)
  disp(sprintf('(dofmrigru1) Could not find any car/ext files in %s',carextDir));
  return
end

% find pdf files
pdfList = getFileList(pdfDir,'pdf');

% get stimfile list
stimfileList = getFileList(stimfileDir,'mat');

% display what we found
dispList(fidList,epiNums,sprintf('Epi scans: %s',fidDir));
dispList(fidList,anatNums,sprintf('Anatomy scans: %s',fidDir));
if senseProcessing
  dispList(fidList,senseRefNums,sprintf('Sense reference scan: %s',fidDir));
  dispList(fidList,senseNoiseNums,sprintf('Sense noise scan: %s',fidDir));
end
dispList(carList,nan,sprintf('Car/Ext files: %s',carextDir));
dispList(pdfList,nan,sprintf('PDF files: %s',pdfDir));
dispList(stimfileList,nan,sprintf('Stimfiles: %s',stimfileDir));

% go find the matching car files for each scan
carMatchNum = getCarMatch(carList,fidList,epiNums);
if isempty(carMatchNum),return,end

% setup directories etc.
doMoveFiles(1,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums);

% make mask dir
if senseProcessing,makeMaskDir(1,fidList,anatNums,senseRefNums);end

% now ask the user if they want to continue, because now we'll actually copy the files and set everything up.
if ~askuser('OK to run above commands?'),return,end

% now do it
fidList = doMoveFiles(0,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums);

% make mask dir
if senseProcessing,makeMaskDir(0,fidList,anatNums,senseRefNums);end

% now run mrInit
disp(sprintf('(dofmrigru1) Setup mrInit for your directory'));
mrInit;

% set up motion comp parameters
for i = 1:numMotionComp
  v = newView;
  [v params] = motionComp(v,[],'justGetParams=1');
  deleteView(v);
  eval(sprintf('save motionCompParams%i params',i));
end

%%%%%%%%%%%%%%%%%%%%%
%    makeMaskDir    %
%%%%%%%%%%%%%%%%%%%%%
function makeMaskDir(justDisplay,fidList,anatNums,senseRefNums)

disp(sprintf('=============================================='));
disp(sprintf('Make empty MLR directory for mask creation    '));
disp(sprintf('=============================================='));

% process ref and anatomy fid to make the mask
command = sprintf('cd Pre;'); myeval(command,justDisplay);% move into Pre
global maskdir;
% location of masks
% make new empty mask directory in Pre
if ~justDisplay
  makeEmptyMLRDir('Mask','description=Empty session for making a mask','defaultParams=1');
else
  disp('makeEmptyMLRDisplay');
end

maskdir = fullfile(pwd,'Mask');

for i = 1:length(senseRefNums)
  disp(sprintf('=============================================='));
  disp(sprintf('Convert ref to nifti and copy into maskdir'));
  disp(sprintf('=============================================='));
  % convert ref to nifti
  if ~justDisplay
    fid2nifti(fidList{senseRefNums(i)}.filename);
  end
  
  % move into emptyMLR directory for mask
  command = sprintf('movefile %s %s;',setext(fidList{senseRefNums(i)}.filename,'hdr'),fullfile(maskdir,'Anatomy'));
  myeval(command,justDisplay);
  command = sprintf('movefile %s %s;',setext(fidList{senseRefNums(i)}.filename,'img'),fullfile(maskdir,'Anatomy'));
  myeval(command,justDisplay);

end  
% move out of Pre
command = sprintf('cd ..;'); myeval(command,justDisplay);

%%%%%%%%%%%%
%% myeval %%
%%%%%%%%%%%%
function myeval(command,justDisplay)

if justDisplay
  disp(command);
else
  eval(command);
end  
%%%%%%%%%%%%%%%%%%%%%
%%   doMoveFiles   %%
%%%%%%%%%%%%%%%%%%%%%
function fidList = doMoveFiles(justDisplay,fidList,carList,pdfList,stimfileList,carMatchNum,epiNums,anatNums,senseNoiseNums,senseRefNums)

% some command names
global epirri;
global postproc;

disp(sprintf('=============================================='));
disp(sprintf('Making directories'));
disp(sprintf('=============================================='));

% list of directories to make
dirList = {'Etc','Pre','Doc','Pre/Aux','Raw','Raw/TSeries','Anatomy','Anal'};

% make them
for i = 1:length(dirList)
  if ~isdir(dirList{i})
    command = sprintf('mkdir(''%s'');',dirList{i});
    if justDisplay,disp(command),else,eval(command);,end
  end
end

% open the logfile
if ~justDisplay
  cd('Pre');
  openLogfile('dofmrigru.log');
  cd('..');
end

disp(sprintf('=============================================='));
disp(sprintf('Copying pdf files'));
disp(sprintf('=============================================='));

% move all the pdf files into the directory
for i = 1:length(pdfList)
  command = sprintf('copyfile %s %s',pdfList{i}.fullfile,fullfile('Doc',pdfList{i}.filename));
  if justDisplay,disp(command),else,eval(command);,end
end

disp(sprintf('=============================================='));
disp(sprintf('Copying stimfiles'));
disp(sprintf('=============================================='));

% move stimfiles
for i = 1:length(stimfileList)
  command = sprintf('copyfile %s %s',stimfileList{i}.fullfile,fullfile('Etc',stimfileList{i}.filename));
  if justDisplay,disp(command),else,eval(command);disp(command);,end
end

if justDisplay,disp(sprintf('+++++++'));end

disp(sprintf('=============================================='));
disp(sprintf('Copying fid files'));
disp(sprintf('=============================================='));

% move fid directories
allUsefulFids = [epiNums anatNums senseNoiseNums senseRefNums];
for i = 1:length(fidList)
  % fix name in our fidList as any . will be replaced by _ in our Pre
  % directory but not in /usr1 or wherever the files were copied from,
  fidList{i}.filename = setext(fixBadChars(stripext(fidList{i}.filename),{'.','_'}),'fid');
  % not a useful scan, put it in Pre/Aux
  if ~any(i==allUsefulFids)
    command = sprintf('copyfile %s %s',fidList{i}.fullfile,fullfile('Pre/Aux',fidList{i}.filename));
    if justDisplay,disp(command);else,eval(command);,end
  %otherwise put it in Pre
  else
    command = sprintf('copyfile %s %s',fidList{i}.fullfile,fullfile('Pre',fidList{i}.filename));
    if justDisplay,disp(command),else,eval(command);,end
  end
end

disp(sprintf('=============================================='));
disp(sprintf('Copying car/ext files'));
disp(sprintf('=============================================='));

% move car/ext files into fid directories
for i = 1:length(carMatchNum)
  command = sprintf('copyfile %s %s',carList{carMatchNum(i)}.fullfile,fullfile('Pre',fidList{epiNums(i)}.filename));
  if justDisplay,disp(command),else,eval(command),end
  command = sprintf('copyfile %s %s',fullfile(carList{carMatchNum(i)}.path,carList{carMatchNum(i)}.extfilename),fullfile('Pre',fidList{epiNums(i)}.filename));
  if justDisplay,disp(command),else,eval(command);,end
end

disp(sprintf('=============================================='));
disp(sprintf('Processing epi files'));
disp(sprintf('=============================================='));

command = sprintf('cd Pre');
if justDisplay,disp(command),else,eval(command),end

% run epirri on all scans and sense noise/ref 
allScanNums = [epiNums senseNoiseNums senseRefNums];
for i = 1:length(allScanNums)
  command = sprintf('status = mysystem(''%s %s 1 1 2 0 1 0'');',epirri,fidList{allScanNums(i)}.filename);
  if justDisplay
    disp(command)
  else
    eval(command);
    if status > 1
        disp(sprintf('(dofmrigru) %s generated error %i when running',epirri,status));
        keyboard
    end    
  end
end

% convert sense noise/ref into edt files
senseNums = [senseNoiseNums senseRefNums];
for i = 1:length(senseNums)
  %make fid to edt
  command = sprintf('mysystem(''%s -intype 0 -outtype 1 %s %s'');',postproc,setext(fidList{senseNums(i)}.filename,'fid'),setext(fidList{senseNums(i)}.filename,'edt'));
  if justDisplay,disp(command),else,eval(command);,end
end

command = sprintf('cd ..');
if justDisplay,disp(command),else,eval(command);,end


disp(sprintf('=============================================='));
disp(sprintf('Making temporary (empty) nifti files for running mrInit'));
disp(sprintf('=============================================='));
for i = 1:length(epiNums)
  srcName = fullfile('Pre',fidList{epiNums(i)}.filename);
  destName = fullfile('Raw/TSeries',setext(fidList{epiNums(i)}.filename,'hdr'));
  if ~justDisplay
    % make a nifti header
    [h info] = fid2niftihdr(srcName);
    % now, grab the data
    if info.accFactor == 1
      d = fid2nifti(srcName);
    % for a sense file we are going to have to grab the reference scan data
    else
      d = fid2nifti(fullfile('Pre',fidList{senseRefNums(1)}.filename));
      % now make the data info a single volume
      d = mean(d(:,:,:,2:end),4);
      % and replicate for the correct number of frames
      d = repmat(d,[1 1 1 h.dim(5)]);
    end
    % write as a nifti file
    cbiWriteNifti(destName,d,h);
  else
    disp(sprintf('Make temporary nifti file for %s in %s',srcName,destName));
  end
end

disp(sprintf('=============================================='));
disp(sprintf('DONE moving files.'));
disp(sprintf('=============================================='));

% convert anatomy scans
for i = 1:length(anatNums)
  disp(sprintf('=========================='));
  disp(sprintf('Convert %s to nifti  ',fidList{anatNums(i)}.filename));
  disp(sprintf('=========================='));
  % convert anatomy to nifti
  command = sprintf('fid2nifti %s;',fullfile('Pre',fidList{anatNums(i)}.filename));
  if justDisplay,disp(command),else,eval(command);end
end


if ~justDisplay
  closeLogfile
end

%%%%%%%%%%%%%%%%%%%%%
%%   getCarmatch   %%
%%%%%%%%%%%%%%%%%%%%%
function carMatchNum = getCarMatch(carList,fidList,epiNums);

carMatchNum = [];

% get the names of the car files
for i = 1:length(carList)
  carNames{i} = carList{i}.filename;
end

% get info from scans
for i = 1:length(epiNums)
  scanNames{i} = stripext(fidList{epiNums(i)}.filename);
  carMatchNames{i} = putOnTopOfList(carNames{min(i,end)},carNames);
  startTime{i} = fidList{epiNums(i)}.startTimeStr;
  endTime{i} = fidList{epiNums(i)}.endTimeStr;
end

% make an mrParamsDialog where you can select
paramsInfo{1} = {'scanNum',1,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',length(epiNums))};
paramsInfo{2} = {'scanName',scanNames,'group=scanNum','type=String','editable=0'};
paramsInfo{3} = {'startTime',startTime,'group=scanNum','type=String','editable=0'};
paramsInfo{4} = {'endTime',endTime,'group=scanNum','type=String','editable=0'};
paramsInfo{5} = {'carFile',carMatchNames,'group=scanNum'};

% put out the dialog
params = mrParamsDialog(paramsInfo,'Match the car file with the scan');
if isempty(params),return;end

% now get the matching numbers
for i = 1:length(params.carFile)
    carMatchNum(i) = find(strcmp(params.carFile(i),carNames));
end


%%%%%%%%%%%%%%%%%%%%%%
%%   dispScanlist   %%
%%%%%%%%%%%%%%%%%%%%%%
function dispStr = dispList(fidList,nums,name)

dispStr = {};
disp(sprintf('============================='));
disp(sprintf('%s',name));
disp(sprintf('============================='));

% if nums is nan, show all files
if isnan(nums),nums = 1:length(fidList);end

for i = 1:length(nums)
  if isfield(fidList{nums(i)},'dispstr')
    disp(fidList{nums(i)}.dispstr);
  else
    disp(fidList{nums(i)}.filename);
  end
end

% empty nums means to display all
if isempty(nums)
  disp('NO MATCHING SCANS');
end

%%%%%%%%%%%%%%%%%%%%
%%   sortFidList  %%
%%%%%%%%%%%%%%%%%%%%
function fidList = sortFidList(fidList)

% sort by time stamp from log file
for i = 1:length(fidList)
  for j = 1:length(fidList)-1
    if fidList{j}.startTime > fidList{j+1}.startTime
      temp = fidList{j};
      fidList{j} = fidList{j+1};
      fidList{j+1} = temp;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%
%%   getCarInfo   %%
%%%%%%%%%%%%%%%%%%%%
function carList = getCarInfo(carList)

disppercent(-inf,'(dofmrigru1) Get car info');
for i = 1:length(carList)
  disppercent(i/length(carList));
  carList{i}.car = readcar(carList{i}.fullfile);
  carList{i}.ext = readext(fullfile(carList{i}.path,carList{i}.extfilename));
  carList{i}.dispstr = sprintf('%s',carList{i}.filename);
end
disppercent(inf);
%%%%%%%%%%%%%%%%%%%%
%%   getFidInfo   %%
%%%%%%%%%%%%%%%%%%%%
function fidList = getFidInfo(fidList)

nVols = [];
disppercent(-inf,'(dofmrigru1) Getting fid info');
for i = 1:length(fidList)
  disppercent(i/length(fidList));
  % read the procpar
  fidList{i}.procpar = readprocpar(fidList{i}.fullfile);
  % get info using fid2xform if we can
  if isfield(fidList{i}.procpar,'ni')
    [fidList{i}.xform fidList{i}.info] = fid2xform(fidList{i}.procpar,-1);
  else
    fidList{i}.xform = [];fidList{i}.info = [];
  end
  % load the logfile
  log = 0;
  if isfile(fullfile(fidList{i}.fullfile,'log'))
    f = fopen(fullfile(fidList{i}.fullfile,'log'));
    thisline = fgets(f);
    [fidList{i}.startTime fidList{i}.startTimeStr] = getDatenumFromLogLine(thisline);
    thisline = fgets(f);
    [fidList{i}.endTime fidList{i}.endTimeStr] = getDatenumFromLogLine(thisline);
    fclose(f);
    log = 1;
  end
  % make the display string
  fidList{i}.dispstr = sprintf('%s: ',fidList{i}.filename);
  fidList{i}.dispstr = sprintf('%s%s->%s',fidList{i}.dispstr,fidList{i}.startTimeStr,fidList{i}.endTimeStr);
  if ~isempty(fidList{i}.info)
    fidList{i}.dispstr = sprintf('%s [%i %i %i %i] [%0.1f %0.1f %0.1f]',fidList{i}.dispstr,fidList{i}.info.dim(1),fidList{i}.info.dim(2),fidList{i}.info.dim(3),fidList{i}.info.dim(4),fidList{i}.info.voxsize(1),fidList{i}.info.voxsize(2),fidList{i}.info.voxsize(3));
     if ~isempty(fidList{i}.info.accFactor)
      fidList{i}.dispstr = sprintf('%s x%i',fidList{i}.dispstr,fidList{i}.info.accFactor);
     end
  end
  fidList{i}.hasCarExt = [0 0];
  fidList{i}.hasPeakFiles = [0 0 0];
  if ~isempty(dir(fullfile(fidList{i}.fullfile,'*.car')))
    fidList{i}.hasCarExt(1) = 1;
  end
  if ~isempty(dir(fullfile(fidList{i}.fullfile,'*.ext')))
    fidList{i}.hasCarExt(2) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'acq.peak.bit'))
    fidList{i}.hasPeakFiles(1) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'cardio.peak.bit'))
    fidList{i}.hasPeakFiles(2) = 1;
  end
  if isfile(fullfile(fidList{i}.fullfile,'respir.peak.bit'))
    fidList{i}.hasPeakFiles(3) = 1;
  end
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getDatenumFromLogLine   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thisDatenum thisDatestr]= getDatenumFromLogLine(thisline)

thisDatenum = 0;thisDatestr = '';
if ~isempty(thisline)
  splitThisLine = strfind(thisline,': ');
  if ~isempty(splitThisLine)
    thisline = thisline(1:splitThisLine-1);
    [dow month day hhmmss year] = strread(thisline,'%s %s %s %s %s');
    thisDatestr = sprintf('%s-%s-%s %s',day{1},month{1},year{1},hhmmss{1});
    thisDatenum = datenum(thisDatestr);
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%%   getSenseNums   %%
%%%%%%%%%%%%%%%%%%%%%%
function [senseNoiseNums senseRefNums] = getSenseNums(fidList)

senseNoiseNums = [];senseRefNums = [];

for i = 1:length(fidList)
  if ~isempty(fidList{i}.info)
    % make sure this is not a 3d acq or a sense
    if ~fidList{i}.info.acq3d && (fidList{i}.info.accFactor == 1)
      % make sure it's not raw
      if fidList{i}.info.dim(3) == length(fidList{i}.info.pss)
        if ~isempty(strfind(fidList{i}.filename,'noise'))
	  senseNoiseNums(end+1) = i;
        elseif ~isempty(strfind(fidList{i}.filename,'ref'))
	  senseRefNums(end+1) = i;
        end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getSenseScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function senseNums = getSenseScanNums(fidList,epiNums)

senseNums = [];
for i = 1:length(epiNums)
  if ~isempty(fidList{epiNums(i)}.info)
    if ~isempty(fidList{epiNums(i)}.info.accFactor)
      senseNums(end+1) = i;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%   getEpiScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function epiNums = getEpiScanNums(fidList,nVolsCutoff)

if ieNotDefined('nVolsCutoff'),nVolsCutoff = 10;end

epiNums = [];
for i = 1:length(fidList)
  if ~isempty(fidList{i}.info)
    if fidList{i}.info.dim(4) > nVolsCutoff
      epiNums(end+1) = i;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getAnatScanNums   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function anatNums = getAnatScanNums(fidList)

anatNums = [];
% look for 3D scans
for i = 1:length(fidList)
  if ~isempty(fidList{i}.info) && fidList{i}.info.processed
    % check for 3d scan
    if fidList{i}.info.acq3d
      % make sure it is not a raw scan
      if fidList{i}.info.dim(3) == length(fidList{i}.info.pss)
        anatNums(end+1) = i;
      end
    else
      % check for an anatomy sounding name (i.e. has the word "anat" in it
      if ~isempty(strfind(lower(fidList{i}.filename),'anat'))
%          keyboard
%          %make sure it is not a raw scan (has been epirri5 processed)
%          command = sprintf('thisanatprocpar = readprocpar(''%s'')',fidList{i}.fullfile);
%          eval(command);
%          if thisanatprocpar.ni ~= 0
            anatNums(end+1) = i;
%         end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%
%%   getFileList   %%
%%%%%%%%%%%%%%%%%%%%%
function [fileList fileListArray] = getFileList(dirname,extList,matchExtList,filenameMatch)

fileList = {};fileListArray = {};

% make into a cell array
extList = cellArray(extList);

if ~ieNotDefined('matchExtList')
  matchExtList = cellArray(matchExtList);
else
  matchExtList = [];
end

% default to no filename matching 
if ieNotDefined('filenameMatch')
  filenameMatch = '';
else
  filenameMatch = cellArray(filenameMatch);
end

% open the directory
dirList = dir(dirname);
if isempty(dirList),return,end

% now go through the directory looking for matches
for i = 1:length(dirList)
  match = 0;
  % skip all . files
  if dirList(i).name(1) == '.',continue,end
  % skip all files that don't match the filename match if specified
  if ~isempty(filenameMatch)
    noFilenameMatch = 0;
    for j = 1:length(filenameMatch)
      if isempty(strfind(lower(dirList(i).name),lower(filenameMatch{j}))),noFilenameMatch=1;,end
    end
    if noFilenameMatch,continue,end
  end
  % get the file extension
  thisExt = getext(dirList(i).name);
  % check for match
   for j = 1:length(extList)
    if strcmp(extList{j},thisExt)
      % found the match
      if ~match
	% keep the name of the file
	fileList{end+1}.filename = dirList(i).name;
	fileListArray{end+1} = dirList(i).name;
	fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	fileList{end}.date = dirList(i).date;
	fileList{end}.datenum = dirList(i).datenum;
	% if we need to find matches
	if ~isempty(matchExtList)
	  fileList{end}.filename = dirList(i).name;
	  fileListArray{end+1} = dirList(i).name;
	  fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	  fileList{end}.path = dirname;
	  fileList{end}.(sprintf('%sfilename',extList{j})) = dirList(i).name;
	  % check for matching file
	  stemName = stripext(fullfile(dirname,dirList(i).name));
	  for k = 1:length(matchExtList)
	    if ~isfile(setext(stemName,matchExtList{j}))
	      disp(sprintf('(dofmrigru1:getFileList) No matching %s file for %s',matchExtList{j},dirList(i).name));
	    else
	      fileList{end}.(sprintf('%sfilename',matchExtList{j})) = setext(dirList(i).name,matchExtList{j});
	    end
	  end
	  
	end
      end
      match = 1;
    end
  end
end

%%%%%%%%%%%%%%%%%%
%%   mysystem   %%
%%%%%%%%%%%%%%%%%%
function status = mysystem(commandName)

% display to buffer
disp('=================================================================');
disp(sprintf('%s',datestr(now)));
disp(commandName);
disp('=================================================================');
% run the command
[status result] = system(commandName);
disp(result);


% write into the logfile
writeLogfile('\n=================================================================\n');
writeLogfile(sprintf('%s\n',datestr(now)));
writeLogfile(sprintf('%s\n',commandName));
writeLogfile('=================================================================\n');
writeLogfile(result);

%%%%%%%%%%%%%%%%%%%%%
%%   openLogfile   %%
%%%%%%%%%%%%%%%%%%%%%
function openLogfile(filename)

global gLogfile;

% open logfile
gLogfile.fid = fopen(filename,'w');
if gLogfile.fid == -1
  disp(sprintf('(dofmrigru) Could not open logfile %s',filename))
  return
end

% remember filename
gLogfile.filename = filename;

%%%%%%%%%%%%%%%%%%%%%%
%%   closeLogfile   %%
%%%%%%%%%%%%%%%%%%%%%%
function closeLogfile

global gLogfile;

% close logfile
if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fclose(gLogfile.fid);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   writeLogfile   %%
%%%%%%%%%%%%%%%%%%%%%%
function writeLogfile(text)

global gLogfile;

if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fprintf(gLogfile.fid,text);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    checkCommands    %
%%%%%%%%%%%%%%%%%%%%%%%
function retval = checkfMRISupportUnitCommands

global epirri;
global postproc;
global sense;

retval = 1;

% commands to check
commandNames = {'epirri','postproc','sense'};
helpFlag = {'','-help',''};
for i = 1:length(commandNames)
  % suse which to tell if we have the command
  [commandStatus commandRetval] = system(sprintf('which %s',eval(commandNames{i})));
  % check for commandStatus error
  if commandStatus~=0
    disp(sprintf('(dofmrigru) Could not find %s command: %s',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmrigru for help setting up your computer'));
    retval = 0;
    return
  end
  % run the command to see what happens
  [commandStatus commandRetval] = system(sprintf('%s %s',eval(commandNames{i}),helpFlag{i}));
  % check for commandStatus error
  if commandStatus>1
    disp(commandRetval);
    disp(sprintf('(dofmrigru) Found %s command: %s, but could not run (possibly missing fink library?)',commandNames{i},eval(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmrigru for help setting up your computer'));
    retval = 0;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%
% isSenseProcessing %
%%%%%%%%%%%%%%%%%%%%%
function senseProcessing = isSenseProcessing(fidList,epiNums)
senseProcessing = 0;

for i = epiNums
  if isfield(fidList{i}.info,'accFactor') && (fidList{i}.info.accFactor>1)
    senseProcessing = 1;
  end
end
