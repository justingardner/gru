classdef GearInvocation < handle
    properties
        gear
        config
        inputs
        destination
        analysisLabel
        tags

        context_
    end
    methods
        function obj = GearInvocation(gear)
            obj.gear = gear;
            obj.config = gear.getDefaultConfig();
            obj.inputs = flywheel.model.JobInputsObject;
            obj.tags = {};
        end
        function obj = updateConfig(obj, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            fields = fieldnames(body);
            for i = 1:numel(fields)
                name = fields{i};
                obj.config.(name) = body.(name);
            end
        end
        function obj = setConfig(obj, key, value)
            obj.config.(key) = value;
        end
        function obj = setAnalysisLabel(obj, label)
            if ~obj.gear.isAnalysisGear()
                throw(MException('ApiClient:apiException', '%s is not an analysis gear!', obj.gear.gear.name));
            end
            obj.analysisLabel = label;
        end
        function obj = addTag(obj, tag)
            obj.tags = [obj.tags, {tag}];
        end
        function obj = addTags(obj, tags)
            obj.tags = [obj.tags, tags];
        end
        function obj = setInput(obj, key, input)
            obj.inputs.(key) = flywheel.Util.toRef(input);
        end
        function obj = setDestination(obj, dest)
            obj.destination = flywheel.Util.toDestination(dest);
        end
        function result = run(obj)
            job = flywheel.model.Job(...
                'gearId', obj.gear.id, ...
                'inputs', obj.inputs, ...
                'destination', obj.destination, ...
                'config', obj.config, ...
                'tags', obj.tags);
            
            if obj.gear.isAnalysisGear()
                % Create a new analysis object using label and job object
                result = obj.addAnalysis(job);
                return
            end

            result = obj.context_.addJob(job);
        end
        function result = proposeBatch(obj, containers, varargin)
            p = inputParser;
            addParameter(p, 'optionalInputPolicy', 'ignored');
            parse(p, varargin{:});

            targets = {};
            for i = 1:numel(containers)
                targets = [targets, struct(flywheel.Util.toRef(containers{i}))];
            end

            if obj.gear.isAnalysisGear()
                analysis = flywheel.model.AnalysisInputAny('label', obj.getAnalysisLabel());
            else
                analysis = [];
            end

            proposal = flywheel.model.BatchProposalInput(...
                'gearId', obj.gear.id, ...
                'config', obj.config, ...
                'tags', obj.tags, ...
                'optionalInputPolicy', p.Results.optionalInputPolicy, ...
                'analysis', analysis, ...
                'targets', targets);

            result = obj.context_.proposeBatch(proposal);
        end
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;
        end
        function result = addAnalysis(obj, job)
            % Ensure a valid analysis label
            label = obj.getAnalysisLabel();

            if isempty(obj.destination)
                throw(MException('ApiClient:apiException', 'Must specify a valid destination to create an analysis!'));
            end

            % Get destination type (from ref)
            destType = [ upper(obj.destination.type(1)) obj.destination.type(2:end) ];
            fn = str2func(sprintf('add%sAnalysis', destType));

            % Ensure destination is unset
            remove(job, 'destination');

            analysis = flywheel.model.AnalysisInput('label', label, 'job', job);
            result = fn(obj.context_, obj.destination.id, analysis);
        end
        function result = getAnalysisLabel(obj)
            if isempty(obj.analysisLabel)
                result = sprintf('%s %s', obj.gear.gear.name, ...
                    datestr(datetime, 'yyyy-mm-dd HH:MM:SS'));
            else
                result = obj.analysisLabel;
            end
        end
    end
end
