classdef Finder < handle
    properties(Hidden)
        context_
        args_
        fn_
    end
    methods
        function obj = Finder(context, method, varargin)
            obj.context_ = context;
            obj.args_ = varargin;
            obj.fn_ = str2func(method);
        end
        function result = subsref(obj, s)
            if numel(s) == 1 && strcmp(s.type, '()')
                % Handle call
                args = horzcat(obj.args_, s.subs);
                result = obj.fn_(obj.context_, args{:});
            else
                result = builtin('subsref', obj, s);
            end
        end
        function result = find(obj, varargin)
            result = obj.find_(false, false, varargin);
        end
        function result = findOne(obj, varargin)
            result = obj.find_(false, true, varargin);
        end
        function result = findFirst(obj, varargin)
            result = obj.find_(true, false, varargin);
        end
        function result = iter(obj, varargin)
            p = inputParser;
            addParameter(p, 'limit', 250);
            parse(p, varargin{:});
            result = obj.iterFind('limit', p.Results.limit);
        end
        function result = iterFind(obj, varargin)
            % Parse out filter strings
            args = obj.makeArgs_(varargin);

            % Ensure that limit is set
            p = inputParser;
            p.StructExpand = false;
            p.KeepUnmatched = true;
            addParameter(p, 'limit', []);
            parse(p, args{:});
            if isempty(p.Results.limit)
                args = horzcat({'limit', 250}, args);
            end
            args = horzcat(obj.args_, args);

            % Generator closure that tracks the last ID
            afterId = '';
            function results = generator()
                currentArgs = horzcat(args, {'afterId', afterId});

                results = obj.fn_(obj.context_, currentArgs{:});
                if ~isempty(results)
                    afterId = results{end}.id;
                end
            end

            % Return a cursor for generator
            result = flywheel.Cursor(@generator);
        end
    end
    methods(Hidden)
        function results = find_(obj, findFirst, findOne, args)
            args = horzcat(obj.args_, obj.makeArgs_(args));
            results = obj.fn_(obj.context_, args{:});

            if findOne
                if numel(results) ~= 1
                    throw(MException('ApiClient:apiException', 'Found %d results instead of 1!', numel(results)));
                end
                results = results{1};
            end
            if findFirst
                if numel(results) > 0
                    results = results{1};
                else
                    results = [];
                end
            end
        end
        function args = makeArgs_(obj, args)
            % Concatenate filter strings with ','
            filter = {};
            operators = {'=', '<', '>'}; 
            for i = 1:numel(args)
                if any(cellfun(@(x) ~isempty(strfind(args{i}, x)), operators))
                    filter = [filter; args{i}];
                else
                    i = i - 1;
                    break
                end
            end
            if ~isempty(filter)
                rest = args(i+1:end);
                filter = strjoin(filter, ',');
                args = horzcat({'filter', filter}, rest);
            end
        end
    end
end
