classdef ResolverOutputMixin < handle
    properties
        context_
    end
    methods(Hidden)
        function setContext_(obj, context)
            obj.context_ = context;

            parent = [];

            % Initialize path files
            path = obj.get('path');
            for i = 1:numel(path)
                node = path{i};
                if isprop(node, 'containerType')
                    if strcmp(node.containerType, 'file')
                        node.parent_ = parent;
                    end
                    parent = node;
                end
            end

            % Initialize child files
            children = obj.get('children');
            for i = 1:numel(children)
                node = children{i};
                if isprop(node, 'containerType') && strcmp(node.containerType, 'file')
                    node.parent_ = parent;
                end
            end
        end
    end
end
