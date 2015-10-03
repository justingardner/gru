% readcar.m
%
%      usage: readcar(filename,<verbose=0>)
%         by: justin gardner
%       date: 07/28/03
%    purpose: Read car file
%
function d = readcar(filename,varargin)

if ~any(nargin==[1 2])
  help readcar;
  return
end

verbose = 0;
getArgs(varargin,{'verbose=0'});

% check which platform we are running (Unix, Intel Mac, Power Mac)
% MAC, MACI, PCWIN, GLNX86
[platform maxsize endian] = computer;

% set filename
d.filename = filename;

% set extension
if ~isfile(filename)
  filename = setext(filename,'car');
end

% open the file
fid = fopen(filename,'r');
if (fid == -1)
  disp(sprintf('ERROR: Could not open %s',filename));
  return
end

% check the first byte to see if we have a car1 or car2 file
if strcmp(endian, 'L')
  [firstByte, n] = fread(fid,1,'uint8=>uint8','ieee-be');
else
  [firstByte, n] = fread(fid,1,'uint8=>uint8');
end

if bitand(firstByte,hex2dec('80')) ~= 0
  if verbose,disp(sprintf('(readcar) Reading Car version 2 file'));end
  d.carVersion = 2;
  d.numChannels = double(bitshift(bitand(firstByte,hex2dec('7C')),-2)+1);
else
  d.carVersion = 1;
  d.numChannels = 16;
  if verbose,disp(sprintf('(readcar) Reading Car version 1 file'));end
end

% seek back to beginning
fseek(fid,0,-1);

% read the file
switch d.carVersion
 case 1
  d = readVersion1(fid,d,endian);
 case 2
  d = readVersion2(fid,d,endian);
end


%%%%%%%%%%%%%%%%%%%%%%
%    readVersion2    %
%%%%%%%%%%%%%%%%%%%%%%
function d = readVersion2(fid,d,endian)

% default fields
d.acq = [];
d.resp = [];
d.cardio = [];
d.channels = [];

% read the data
if strcmp(endian, 'L')
    [raw, n] = fread(fid,inf,'uint16=>uint16','ieee-be');
else
    [raw, n] = fread(fid,inf,'uint16=>uint16');
end

% the n should be an even multiple of number of channels
if rem(n,d.numChannels) ~= 0
  disp(sprintf('(readcar) **** Car file length is incorrect. %i samples does not divide evenly by number of channels: %i ****',n,d.numChannels));
  return
end

% get number of samples
numSamples = n/d.numChannels;

% reshape date to get each signal
raw = reshape(raw,d.numChannels,numSamples);

% convert the acquisition channel.
% lower 10 bits are the 10 1 ms samples.
numAcqSamples = 10;
bitMask = bitshift(1,numAcqSamples)-1;
% mask out the data bits and turn into a string of 0/1 for each sample
% fliplr since the lowest bit corresponds to the earliest time
acq = fliplr(dec2bin(bitand(raw(1,:),bitMask),numAcqSamples));
% convert dimensions 
d.acq = double(reshape(acq',1,numAcqSamples*numSamples)=='1');

% now, reread the data this time as signed integers
fseek(fid,0,-1);
if strcmp(endian, 'L')
    [raw, n2] = fread(fid,inf,'int16=>int16','ieee-be');
else
    [raw, n2] = fread(fid,inf,'int16=>int16');
end

if n ~= n2
  disp(sprintf('(readcar) **** Error reading file, should have read %i ints, but read %i instaed ****',n,n2));
  return
end

% reshape date to get each signal
raw = reshape(raw,d.numChannels,numSamples);

% now get cardio and resp, these are stored as int16 in rows 2 + 3
d.cardio = double(raw(2,:));
d.resp = double(raw(3,:));

% and get the rest of the channels
d.channels = double(raw(4:end,:));

%close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%
%    readVersion1    %
%%%%%%%%%%%%%%%%%%%%%%
function d = readVersion1(fid,d,endian)

% read the data
if strcmp(endian, 'L')
    [raw, n] = fread(fid,inf,'uint32=>uint32','ieee-be');
else
    [raw, n] = fread(fid,inf,'uint32=>uint32');
end

% close the file
fclose(fid);

% get data points
d.acq = double(bitand(bitshift(raw,-28),hex2dec('F')));
d.cardio = bitand(bitshift(raw,-16),hex2dec('FFF'));
d.resp = bitand(raw,hex2dec('FFF'));

% for each value greater than 0x7FF replace with -0xFFF+value
d.resp = double(d.resp>hex2dec('7FF')).*-double(hex2dec('FFF'))+double(d.resp);
d.cardio = double(d.cardio>hex2dec('7FF')).*-double(hex2dec('FFF'))+double(d.cardio);
%scale by 5/0x7FF
d.resp = 5*d.resp/double(hex2dec('7ff'));
d.cardio = 5*d.cardio/double(hex2dec('7ff'));

% make into row vector
d.acq = d.acq';
d.cardio = d.cardio';
d.resp = d.resp';
