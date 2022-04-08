classdef AnalysisMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.InfoMethods ...
        & flywheel.mixins.NoteMethods ...
        & flywheel.mixins.FileMethods ...
        & flywheel.mixins.DownloadMethods
    properties
        containerType_ = 'Analysis';
        fileGroup_ = 'Output';
    end
    methods
        function [returnData, resp] = uploadFile(obj, file)
            [returnData, resp] = obj.invokeContainerApi('uploadOutputToAnalysis', obj.get('id'), file);
        end
        function [returnData, resp] = uploadOutput(obj, file)
            [returnData, resp] = obj.invokeContainerApi('uploadOutputToAnalysis', obj.get('id'), file);
        end
    end
end