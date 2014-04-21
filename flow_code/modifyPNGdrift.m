close all
clear all

load('pngNames.mat')

noPic = length(pngNames);

picDim = [276,511];

verticalDrift = 36.1507; %5.3;
horizontalDrift = 34.7595; %4;

%mean(mean(uv))(1) = 34.7595, (2) = -36.1507
%in plotflow, meshgrid x dimension is size(uv)(2) and y dimension is
%size(uv)(1)

%rescale average drift found from OFA to match im1 and im2. Then shift im2
%accordingly and plot on top of each other to check method.

%then 

%then get average drift for each pair of images, and use this

drift = round([linspace(0,verticalDrift*noPic,noPic);linspace(0,horizontalDrift*noPic,noPic)]);

blackPicture = uint8(zeros(picDim(1)+max(drift(1,:)) , picDim(2)+max(drift(2,:)),3));

for pic = 1:noPic
    pngBefore = imread(pngNames(pic,:));
    
    pngAfter = blackPicture;
    pngAfter(drift(1,end)-drift(1,pic)+1:drift(1,end)-drift(1,pic)+picDim(1) , drift(2,end)-drift(2,pic)+1:drift(2,end)-drift(2,pic)+picDim(2),:) = pngBefore;
    
    imwrite(pngAfter, ['new_' pngNames(pic,:)])
end