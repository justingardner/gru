function mlrMRIqc(group,subj)

% see: https://mriqc.readthedocs.io/en/stable/docker.html#docker
% standard code line:
%   docker run -it --rm -v <bids_dir>:/data:ro -v <output_dir>:/out poldracklab/mriqc:latest /data /out participant --participant_label 001 002 003

folder = fullfile('~/data',group,subj,'bids/');
fout = fullfile(folder,'mriqc');

command = sprintf('docker run -it --rm -v %s:/data:ro -v %s:/out poldracklab/mriqc:latest /data /out participant --participant_label 01',folder,fout);

system(command);