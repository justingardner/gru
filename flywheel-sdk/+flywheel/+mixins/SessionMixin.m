classdef SessionMixin < flywheel.mixins.ContainerBase ...
        & flywheel.mixins.TagMethods ...
        & flywheel.mixins.InfoMethods ...
        & flywheel.mixins.NoteMethods ...
        & flywheel.mixins.FileMethods ...
        & flywheel.mixins.TimestampMethods ...
        & flywheel.mixins.DownloadMethods ...
        & flywheel.mixins.AnalysisMethods
    properties
        containerType_ = 'Session';
        acquisitions_ = false;
        analyses_ = false;
    end
    properties(Dependent)
        acquisitions
        ageYears
        ageMonths
        ageWeeks
        ageDays
    end
    methods
        function acquisitions = get.acquisitions(obj)
            acquisitions = obj.getChildren('Acquisitions');
        end
        function analyses = get_analyses(obj)
            analyses = obj.getChildren('Analyses');
        end
        function [returnData, resp] = addAcquisition(obj, varargin)
            [returnData, resp] = obj.addChild('Acquisition', varargin{:});
        end
        function result = get.ageYears(obj)
            result = obj.get('age');
            if ~isempty(result)
                result = flywheel.Util.secondsToYears(result);
            end
        end
        function result = get.ageMonths(obj)
            result = obj.get('age');
            if ~isempty(result)
                result = flywheel.Util.secondsToMonths(result);
            end
        end
        function result = get.ageWeeks(obj)
            result = obj.get('age');
            if ~isempty(result)
                result = flywheel.Util.secondsToWeeks(result);
            end
        end
        function result = get.ageDays(obj)
            result = obj.get('age');
            if ~isempty(result)
                result = flywheel.Util.secondsToDays(result);
            end
        end
    end
end