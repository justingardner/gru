

% author: steeve laquitaine
%   date: 151202
%purpose: simulate random values instances for testing classification
%
%  usage:
%
%           c = simInstances({'simV1','simMT'},2,[20 30],[10 20])
%           c = leaveOneOut(c);

function c = simInstances(rois,nClasses,nvoxs,nins)

%# of rois
nRois = length(nvoxs);

%make instances for each roi and 
%class
for roi = 1 : nRois
    %name roi
    c{roi}.name = rois{roi};
    for class = 1 : nClasses        
        c{roi}.classify.instances{class} = rand(nins(class),nvoxs(roi));
    end
end