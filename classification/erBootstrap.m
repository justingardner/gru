function d = erBootstrap(varargin)
% d = erBootstrap(d)
%
% Adds a field d.p which contains the bootstrapped p values for each voxel
% in each of the loaded ROIs in d (use loadroi first). 

%% Get defaults
stimVol = [];
tSeries = [];
iterations = 0;
time = 0;
d = [];
stimVolString='';
view=[];
getArgs('view=[]','d=[]','iterations=-1','time=60','stimVolString=[]');

%% Setup & checks

if isempty(view) || isempty(d)
    help erBootstrap;
    return
end

time0 = -1;
if iterations>-1
    % use iterations
    time = 0;
else
    time0 = mglGetSecs;
    iterations = 10;
end

if isempty(stimVolString)
    disp('(erBootstrap) Defaulting to deconvolution of all trials. Did you mean to set stimVolString?');
    stimVolString = '_all_';
end
[stimVol, stimNames, ~] = getStimvol(view,stimVolString);

% original scm
scm = makescm(view,d.hdrlen,0,stimVol);

%% get the stimVol parameters
diffs = diff(stimVol{1});
nsv = length(stimVol{1});

%% setup parallel pool
po = gcp('nocreate');
if isempty(po)
    po = parpool('local',6);
end

if isempty(po)
    disp('(erBootstrap) Cluster failed to initialize, this code could take a LONG time');
end


%% run bootstrap
tr2 = cell(1,length(d.roi));
while true % run forever, until either we run out of iterations or we run out of time
    % initialize r2
    r2 = cell(1,length(d.roi));
    for ri=1:length(d.roi), end
    % iterate
        disppercent(-1/length(d.roi));
    for ri = 1:length(d.roi)
        r = d.roi{ri};
        
        % compute true r2
        decon = getr2timecourse(r.tSeries,d.nhdr,d.hdrlen,scm,d.tr,0);
        tr2{ri} = decon.r2;
        clear decon

        r2s=zeros(size(d.roi{ri}.tSeries,1),iterations);
        disp(sprintf('Calculating R^2 for ROI: %i\n',ri));
        decon = {}; i_stimVol = {};
        
        td = [];
        td.concatInfo = viewGet(view,'concatInfo');
        td.dim = viewGet(view,'scanDims');
        td.hdrlen = d.hdrlen;
        
        tdo = {};
        tic
        parfor i = 1:iterations
            % generate scm for this iteration
            i_stimVol{i} = cumsum([1 diffs(randperm(length(diffs)))]);
            tdo{i} = makescm(td,td.hdrlen,0,{i_stimVol{i}});
            i_scm = tdo{i}.scm;
            tdo{i} = struct;
            i_stimVol{i} = struct;

            decon{i} = getr2timecourse(r.tSeries,d.nhdr,d.hdrlen,i_scm,d.tr,0);
            r2s(:,i) = decon{i}.r2;
%             if mod(i,10)==0,disppercent(i/iterations);end
            decon{i} = struct;
            if mod(i,100)==0, disp(sprintf('Finished block: %i/%i',i,iterations));end
        end
        toc
        
        r2{ri} = r2s;
        disppercent(ri/length(d.roi));
    end
    % check while loop constraints
    if time0 == -1 || mglGetSecs > time0 + time*60
        break
    end
end

%% compute p-value of this r2

% find p values
p = {};
for ri = 1:length(d.roi)
    rr = sort(r2{ri},2);
    tr = tr2{ri};
    
    ps = zeros(size(tr));
    for i=1:size(rr,1)
        perc = find(tr(i)<rr(i,:),1)/iterations;
        if isempty(perc)
            ps(i,1) = 0;
        else
            ps(i,1) = 1-perc;
        end
    end
    p{ri}=ps;
end

d.p = p;

%% compute p-value range for each roi
disp('(erBootstrap) p-value range skipped');

%%
delete(po);