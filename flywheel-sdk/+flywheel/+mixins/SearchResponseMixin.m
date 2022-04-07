classdef SearchResponseMixin < handle
    properties
        context_
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;

            % Also initialize file parent
            file = obj.get('file');
            if ~isempty(file)
                file.parent_ = obj.get('parent');
            end
        end
    end
end
