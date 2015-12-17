

% author: steeve laquitaine
%   date: 151202
%purpose: simulate instances to test classification
%
%  usage:
%
% ex: 1
%           c=simInstances({'simV1','simMT'},2,[20 30],[10 20])
%           c=leaveOneOut(c,'permutation=1');
%
%
% ex:2
%           %simulate two 2D-Gaussian-clusters of instances
%           c=simInstances({'simV1','simMT'},2,[],[50 50],'type=a2Dclusters')                   
%           c = c{1}.classify.instances;
%           c=leaveOneOut(c,'permutationBal=1');
%           %plot
%           classs={'r','b'};
%           for i=1:length(c{1}.classify.instances)
%              hold on
%              plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%              '.','color',classs{i})
%           end
% 
% ex:3
%           %2D-Gaussian-clusters 
%           c=simInstances({'simV1','simMT'},2,[],[10 15],'type=a2Dclusters')         
%           c=leaveOneOut(c,'balancByBootSt=1');           
%           %plot
%           classs={'r','b'};
%           nClass = length(c{1}.classify.instances);
%           for i=1:nClass
%              hold on
%              plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%              '.','color',classs{i})
%           end
%
% ex4:
%           %2D-Gaussian-clusters 
%           c=simInstances({'simV1','simMT'},2,[],[10 15],'type=a2Dclusters')
%           %classify
%           c=leaveOneOut(c,'balancByRemovI=1');
%           %plot
%           classs={'r','b'};
%           nClass=length(c{1}.classify.instances);
%           for i=1:nClass
%               hold on
%               plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%               '.','color',classs{i})
%           end



function c = simInstances(rois,nClasses,nvoxs,nins,varargin)

type=[];
getArgs(varargin,{'type=rnd'});

%# of rois
nRois = length(rois);

%(default) case random instance values
if strcmp(type,'rnd')
    %make instances for each roi and
    %class
    for roi = 1 : nRois
        %name roi
        c{roi}.name = rois{roi};
        for class = 1 : nClasses
            c{roi}.classify.instances{class} = rand(nins(class),nvoxs(roi));
        end
    end
end

%case 2D cluster
%2 voxels
mu = [1 -1];
sigma = [.9 .4; .4 .3];
if strcmp(type,'a2Dclusters')
    %make instances for each roi and
    %class
    for roi = 1 : nRois
        %name roi
        c{roi}.name = rois{roi};
        for class = 1 : nClasses
            m = mu + (class - 1)*10;     
            c{roi}.classify.instances{class} = mvnrnd(m,sigma,nins(class))*100;
        end
    end
end