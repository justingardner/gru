classdef TagMethods < handle
    methods
        function [returnData, resp] = addTag(obj, tag)
            [returnData, resp] = obj.invokeContainerApi('add%sTag', obj.get('id'), tag);
        end
        function [returnData, resp] = renameTag(obj, tag, newTag)
            [returnData, resp] = obj.invokeContainerApi('rename%sTag', obj.get('id'), tag, newTag);
        end
        function [returnData, resp] = deleteTag(obj, tag)
            [returnData, resp] = obj.invokeContainerApi('delete%sTag', obj.get('id'), tag);
        end
    end
end
