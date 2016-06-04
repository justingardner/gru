function [h,m] = csvreadh( filename, delim )
%CSVREADH Read a comma separated value file with header.
%   [H,M] = CSVREADH('FILENAME') reads a comma separated value formatted file
%   FILENAME.  The result data is returned in M, the header in H. 
%   The file can only contain numeric values as data and a string for 
%   the header.

% Validate input args
if nargin==0
    error(nargchk(1,1,nargin,'struct')); 
end

% Get Filename
if ~ischar(filename)
    error('csvreadh:FileNameMustBeString', ...
        'Filename must be a string.'); 
end

% Make sure file exists
if exist(filename,'file') ~= 2 
    error('csvreadh:FileNotFound',...
    'File not found.');
end

if nargin==1
    delim = ',';
end

% open input file
file = fopen( filename );
line = fgetl( file );
h = regexp( line, delim, 'split' );

m = [];
% dan code
% we have the header, so instead of creating a GIANT array one flip at a
% time, let's initialize it to have 100,000 empty slots and let matlab fill
% it. If it runs out of space we'll add another 100,000. At the end we will
% just chop off the extra.
% Much faster than the stupid code the person who wrote this put here.
m = zeros(100000,size(h,2));
maxcount = 100000;
count = 1;
while 1
    line = fgetl( file );
    if ~ischar(line), break, end
    if count>maxcount
        maxcount = maxcount+100000;
        m = [m; zeros(100000,size(h,2))];
    end
    m(count,:) = str2double(regexp( line, ',', 'split' ));
    count = count+1;
end

m = m(1:count,:);

fclose(file);