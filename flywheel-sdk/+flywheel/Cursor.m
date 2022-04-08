classdef Cursor < handle
    properties
        generator_
        pending_
        index_
    end
    methods
        function obj = Cursor(generator)
            obj.generator_ = generator;
            obj.pending_ = true;
            obj.index_ = 1;
        end
        function result = hasNext(obj)
            % Refill guarantees that pending_ will be an array, or false
            obj.refill();
            result = ~islogical(obj.pending_);
        end
        function result = next(obj)
            % Refill guarantees that pending_ will be an array, or false
            obj.refill();
            if islogical(obj.pending_)
                result = [];
            else
                result = obj.pending_{obj.index_};
                obj.index_ = obj.index_ + 1;
            end
        end
    end
    methods(Hidden)
        function refill(obj)
            % Fetch if obj.pending_ is true, or we've consumed every element in pending_
            fetch = false;
            if islogical(obj.pending_)
                fetch = obj.pending_;
            elseif obj.index_ >= numel(obj.pending_)
                fetch = true;
            end

            if fetch
                obj.index_ = 1;
                obj.pending_ = obj.generator_();
                
                % If we get an empty result, there are no more elements
                if isempty(obj.pending_)
                    obj.pending_ = false;
                end
            end
        end
    end
end