% readcar.m
%
%      usage: readcar.m()
%         by: justin gardner
%       date: 07/28/03
%    purpose: Read car file
%
function d = readcar(filename)

if (nargin ~= 1)
  help readcar;
  return
end

% check which platform we are running (Unix, Intel Mac, Power Mac)
% MAC, MACI, PCWIN, GLNX86
[platform maxsize endian] = computer;

% set filename
d.filename = filename;

% open the file
fid = fopen(filename,'r');
if (fid == -1)
  disp(sprintf('ERROR: Could not open %s',filename));
  return
end

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
