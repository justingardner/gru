%%
cdir = '/Users/dan/data/cogneuro/s032420160121/Etc';
files = dir(cdir);

files = files(3:end);

stims = {};
for fi = 1:length(files)
    stims{end+1} = pRFGetStimImageFromStimfile(fullfile(cdir,files(fi).name));
end

files = files(1:6);
%% Compare mod vs original

for i = 1:3
    im1 = stims{(i-1)*2+1}.im;
    im2 = stims{(i-1)*2+2}.im;
    imdiff = im1-im2;
    error = sum(sum(sum(imdiff,3),2));
    disp(error);
end

%% Compare pair-wise
files = files([1 3 5]);
stims = stims([1 3 5]);

%% Comparison

pairs = [1 2
 1 3
 2 3];

for p = 1:3
    pair = pairs(p,:);
    im1 = stims{pair(1)}.im;
    im2 = stims{pair(2)}.im;
    imdiff = im1-im2;
    error = sum(sum(sum(imdiff,3),2));
    disp(error);
end

%% Vis difference
figure
time = size(im1,3);
for p = 1:3
    pair = pairs(p,:);
    im1 = stims{pair(1)}.im;
    im2 = stims{pair(2)}.im;
    imdiff = im1-im2;
    for t = 1:time
        subplot(131);
        imshow(squeeze(im1(:,:,t)));
        subplot(132);
        imshow(squeeze(im2(:,:,t)));
        subplot(133);
        imshow(squeeze(imdiff(:,:,t)));
        pause(.01);
    end
end