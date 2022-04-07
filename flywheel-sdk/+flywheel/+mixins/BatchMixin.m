classdef BatchMixin < handle
    properties
        context_
    end
    methods
        function result = cancel(obj)
            result = obj.context_.cancelBatch(obj.id);
        end
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;
        end
    end
end