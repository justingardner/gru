classdef GearMixin < handle
    properties
        context_
    end
    methods
        function result = createInvocation(obj)
            result = flywheel.mixins.GearInvocation(obj);
            result.setContext_(obj.context_);
        end
        function result = isAnalysisGear(obj)
            result = strcmp(obj.category, 'analysis');
        end
        function obj = printDetails(obj)
            gear = obj.gear;

            % Name + Description
            if ~isempty(gear.label)
                fprintf('%s\n', gear.label);
            else
                fprintf('%s\n', gear.name);
            end

            if ~isempty(gear.description)
                fprintf('\n%s\n', gear.description)
            end

            % Attributes
            fprintf(    'Name:       %s\n', gear.name);
            fprintf(    'Version:    %s\n', gear.version);
            fprintf(    'Category:   %s\n', obj.category);

            if ~isempty(gear.author)
                fprintf('Author:     %s\n', gear.author);
            end
            if ~isempty(gear.maintainer)
                fprintf('Maintainer: %s\n', gear.maintainer);
            end
            if ~isempty(gear.url)
                fprintf('URL:        %s\n', gear.url);
            end
            if ~isempty(gear.source)
                fprintf('Source:     %s\n', gear.source);
            end
            fprintf('\n');

            % Inputs
            inputs = struct(gear.inputs);
            if ~isempty(inputs)
                fprintf('Inputs:\n');

                fields = fieldnames(inputs);
                for i = 1:numel(fields)
                    name = fields{i};
                    field = inputs.(name);

                    if islogical(field.optional) && field.optional == false
                        opt = 'optional';
                    else
                        opt = 'required';
                    end

                    fprintf('  %s (%s, %s)\n', name, field.base, opt);
                    if ~isempty(field.description)
                        fprintf('    %s\n', field.description);
                    end
                end
                fprintf('\n');
            end

            % Configuration values
            configKeys = gear.config.keys();
            if ~isempty(configKeys)
                fprintf('Configuration:\n');

                for i = 1:numel(configKeys)
                    name = configKeys{i};
                    field = gear.config.(name);

                    if ~isempty(field.default)
                        dfltStr = strip(evalc('disp(field.default)'));
                        dflt = sprintf(', default: %s', dfltStr);
                    else
                        dflt = '';
                    end

                    fprintf('  %s (%s%s)\n', name, field.type, dflt);
                    if ~isempty(field.description)
                        fprintf('    %s\n', field.description);
                    end
                end
            end
        end
        function result = getDefaultConfig(obj)
            config = struct(obj.gear.config);
            result = flywheel.model.JobConfig;
            fields = fieldnames(config);
            for i = 1:numel(fields)
                name = fields{i};
                field = config.(name);
                result.(name) = field.default;
            end
        end
        function result = run(obj, varargin)
            p = inputParser;
            addParameter(p, 'config', []);
            addParameter(p, 'analysisLabel', []);
            addParameter(p, 'tags', []);
            addParameter(p, 'destination', []);
            addParameter(p, 'inputs', []);
            p.KeepUnmatched = true;

            parse(p, varargin{:});

            invocation = obj.createInvocation();

            if ~isempty(p.Results.config)
                flywheel.Util.applyFnToStructOrCells(@(key, value) invocation.setConfig(key, value), p.Results.config);
            end
            if ~isempty(p.Results.analysisLabel)
                invocation.setAnalysisLabel(p.Results.analysisLabel);
            end
            if ~isempty(p.Results.tags)
                invocation.addTags(p.Results.tags);
            end
            if ~isempty(p.Results.destination)
                invocation.setDestination(p.Results.destination);
            end

            % Set inputs
            if ~isempty(p.Results.inputs)
                flywheel.Util.applyFnToStructOrCells(@(key, value) invocation.setInput(key, value), p.Results.inputs);
            end

            % Convenience to set further inputs from unmatched values
            if ~isempty(p.Unmatched)
                fields = fieldnames(p.Unmatched);
                for i = 1:numel(fields)
                    name = fields{i};
                    invocation.setInput(name, p.Unmatched.(name));
                end
            end

            result = invocation.run();
        end
        function result = proposeBatch(obj, varargin)
            p = inputParser;
            addRequired(p, 'containers');
            addParameter(p, 'config', []);
            addParameter(p, 'analysisLabel', []);
            addParameter(p, 'tags', []);
            addParameter(p, 'optionalInputPolicy', 'ignored');

            parse(p, varargin{:});

            invocation = obj.createInvocation();

            if ~isempty(p.Results.config)
                flywheel.Util.applyFnToStructOrCells(@(key, value) invocation.setConfig(key, value), p.Results.config);
            end
            if ~isempty(p.Results.analysisLabel)
                invocation.setAnalysisLabel(p.Results.analysisLabel);
            end
            if ~isempty(p.Results.tags)
                invocation.addTags(p.Results.tags);
            end

            result = invocation.proposeBatch(p.Results.containers, 'optionalInputPolicy', p.Results.optionalInputPolicy);
        end
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;
        end
    end
end
