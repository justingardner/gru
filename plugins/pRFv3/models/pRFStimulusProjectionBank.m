classdef pRFStimulusProjectionBank < handle
% pRFStimulusProjectionBank
%
% Immutable unique-frame storage shared by prepared stimulus runs. Handle
% identity gives pRFProjectStimulus a collision-proof, constant-time cache
% key, including across repeated preparation calls and function clears.

properties (SetAccess = private)
  Images
  NumBytes
end

methods
  function bank = pRFStimulusProjectionBank(images)
    bank.Images = images;
    imageInfo = whos('images');
    bank.NumBytes = double(imageInfo.bytes);
  end
end
end
