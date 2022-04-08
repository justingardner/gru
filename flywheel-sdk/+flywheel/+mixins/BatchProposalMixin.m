classdef BatchProposalMixin < handle
    properties
        context_
    end
    methods
        function result = run(obj)
            result = obj.context_.startBatch(obj.id);
        end
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;
        end
    end
end