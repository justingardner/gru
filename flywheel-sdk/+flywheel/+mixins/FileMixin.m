classdef FileMixin < flywheel.mixins.ContainerBase
    properties(Hidden)
        containerType_ = 'File';
        parent_
    end
    properties(Dependent)
        url
        parent
    end
    methods
        function result = get.url(obj)
            result = obj.parent_.getFileDownloadUrl(obj.get('name'));
        end
        function result = get.parent(obj)
            result = obj.parent_;
        end
        function resp = download(obj, destFile, varargin)
            resp = obj.parent_.downloadFile(obj.get('name'), destFile, varargin{:});
        end
        function resp = read(obj, varargin)
            resp = obj.parent_.readFile(obj.get('name'), varargin{:});
        end
        function resp = update(obj, varargin)
            resp = obj.parent_.updateFile(obj.get('name'), varargin{:});
        end
        function resp = replaceInfo(obj, info)
            resp = obj.parent_.replaceFileInfo(obj.get('name'), info);
        end
        function resp = updateInfo(obj, varargin)
            resp = obj.parent_.updateFileInfo(obj.get('name'), varargin{:});
        end
        function resp = deleteInfo(obj, varargin)
            resp = obj.parent_.deleteFileInfo(obj.get('name'), varargin{:});
        end
        function resp = replaceClassification(obj, varargin)
            resp = obj.parent_.replaceFileClassification(obj.get('name'), varargin{:});
        end
        function resp = updateClassification(obj, classification)
            resp = obj.parent_.updateFileClassification(obj.get('name'), classification);
        end
        function resp = deleteClassification(obj, classification)
            resp = obj.parent_.deleteFileClassification(obj.get('name'), classification);
        end
        function result = ref(obj)
            result = flywheel.model.FileReference(...
                'type', obj.parent_.containerType, ...
                'id', obj.parent_.id, ...
                'name', obj.name);
        end
        function resp = getZipInfo(obj)
            resp = obj.parent_.getFileZipInfo(obj.get('name'));
        end
        function resp = downloadZipMember(obj, destFile, memberPath, varargin)
            resp = obj.parent_.downloadFileZipMember(obj.get('name'), memberPath, destFile, varargin{:});
        end
        function resp = readZipMember(obj, memberPath, varargin)
            resp = obj.parent_.readFileZipMember(obj.get('name'), memberPath, varargin{:});
        end
    end
end