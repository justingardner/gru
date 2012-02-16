% epibsi.m
%
%        $Id:$ 
%      usage: epibsi()
%         by: justin gardner
%       date: 02/16/12
%    purpose: 
%
function retval = epibsi(filename)

% check arguments
if ~any(nargin == [1])
  help epibsi
  return
end

currentPath = pwd;

filename = setext(filename,'fid');
[x i] = fid2xform(filename);
if isempty(i),return,end

% check that it is an epi
if ~i.isepi
  disp(sprintf('(epibsi) Fid is not an epi. Not running epibsi'));
  return
end

% see if it is compressed
if ~i.compressedFid
  disp(sprintf('(epibsi) Fid is not compressed. No need to run epibsi'));
  return
end

% ok make a filename for the copy
tempFilename = sprintf('~/Desktop/%s_deleteme.fid',stripext(filename));
tempFilename = mlrReplaceTilde(tempFilename);
if isfile(tempFilename)
  disp(sprintf('(epibsi) Filename %s already exists, aborting',tempFilename));
  return
end

% copy the file
disp(sprintf('(epibsi) Copying to temporary folder %s',tempFilename));
system(sprintf('cp -r %s %s',filename,tempFilename));

% run epibsi on the temporary folder
system(sprintf('epibsi %s',tempFilename));

% and display it
mlrVol(tempFilename);















