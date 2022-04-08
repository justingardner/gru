classdef DownloadMethods < handle
    methods
        function summary = downloadTar(obj, varargin)
            % Download the container as a tarball to dest_file.
            %
            % Parameters:
            %   destFile: (required) The destination file on disk
            %   includeTypes: The optional list of types to include in the download (e.g. ['nifti'])
            %   excludeTypes: The optional list of types to exclude from the download (e.g. ['dicom'])
            summary = obj.context_.downloadTar(obj, varargin{:});
        end
    end
end
