classdef ContainerBase < handle
    properties
        context_
    end
    properties(Dependent)
        containerType
        localCreated
        localModified
    end
    methods
        function result = get.containerType(obj)
            result = lower(obj.containerType_);
        end
        function [returnData, resp] = update(obj, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            [returnData, resp] = obj.invokeContainerApi('modify%s', obj.get('id'), body);
        end
        function [returnData, resp] = reload(obj)
            [returnData, resp] = obj.invokeContainerApi('get%s', obj.get('id'));
        end
        function result = ref(obj)
            result = flywheel.model.ContainerReference('type', obj.containerType, 'id', obj.id);
        end
        function result = get.localCreated(obj)
            result = obj.localizeDate('created');
        end
        function result = get.localModified(obj)
            result = obj.localizeDate('modified');
        end
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;
        end
        function [returnData, resp] = invokeContainerApi(obj, name, varargin)
            fn = str2func(sprintf(name, obj.containerType_));
            [returnData, resp] = fn(obj.context_, varargin{:});
        end
        function result = getChildren(obj, childName)
            varName = strcat(lower(childName), '_');
            if islogical(obj.(varName)) && obj.(varName) == false
                fname = sprintf('get%s%s', obj.containerType_, childName);
                obj.(varName) = flywheel.Finder(obj.context_, fname, obj.get('id'));
            end
            result = obj.(varName);
        end
        function [returnData, resp] = addChild(obj, childName, varargin)
            body = flywheel.mixins.ContainerBase.structFromArgs(varargin);
            body.(obj.containerType) = obj.id;
            fn = str2func(sprintf('add%s', childName));
            [childId, resp] = fn(obj.context_, body);
            if ~isempty(childId)
                getFn = str2func(sprintf('get%s', childName));
                returnData = getFn(obj.context_, childId);
            else
                returnData = [];
            end
        end
        function result = localizeDate(obj, key)
            result = obj.get(key);
            if ~isempty(result)
                result = datetime(result, 'TimeZone', 'local');
            end
        end
    end
    methods(Static, Hidden)
        function body = structFromArgs(varargs)
            if numel(varargs) == 1
                body = varargs{1};
            else
                p = inputParser;
                p.StructExpand = false;
                p.KeepUnmatched = true;
                parse(p, varargs{:});
                body = p.Unmatched;
            end

            if isempty(body)
                throw(MException('ApiClient:inputError', 'Must provide a body!'));
            end
        end
    end
end