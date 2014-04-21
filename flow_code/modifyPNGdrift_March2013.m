close all
clear all

load('pngNames.mat')

noPic = length(pngNames);

picDim = [276,511];

%verticalDrift = 36.1507; %5.3;
%horizontalDrift = 34.7595; %4;

%mean(mean(uv))(1) = 34.7595, (2) = -36.1507
%in plotflow, meshgrid x dimension is size(uv)(2) and y dimension is
%size(uv)(1)

%rescale average drift found from OFA to match im1 and im2. Then shift im2
%accordingly and plot on top of each other to check method.

%then 

%then get average drift for each pair of images, and use this

drift = round([0 cumsum(average_drift(:,2))' ; 0 cumsum(average_drift(:,1))']);
%round([linspace(0,verticalDrift*noPic,noPic);linspace(0,horizontalDrift*noPic,noPic)]);

blackPicture = uint8(zeros(picDim(1)+max(drift(1,:)) , picDim(2)+max(drift(2,:)),3));

for pic = 1:noPic-1
    pngBefore = imread(pngNames(pic,:));
    
    pngAfter = blackPicture;
    pngAfter(drift(1,end)-drift(1,pic)+1:drift(1,end)-drift(1,pic)+picDim(1) , drift(2,end)-drift(2,pic)+1:drift(2,end)-drift(2,pic)+picDim(2),:) = pngBefore;
    
    imwrite(pngAfter, ['new_' pngNames(pic,:)])
end



%%%%%%%%%%%%%%%%%%%%%%Removing the full deformation
close all
clear all

load('pngNames.mat')

noPic = length(pngNames);

picDim = [276,511];

%verticalDrift = 36.1507; %5.3;
%horizontalDrift = 34.7595; %4;

%mean(mean(uv))(1) = 34.7595, (2) = -36.1507
%in plotflow, meshgrid x dimension is size(uv)(2) and y dimension is
%size(uv)(1)

%rescale average drift found from OFA to match im1 and im2. Then shift im2
%accordingly and plot on top of each other to check method.

%then 

%then get average drift for each pair of images, and use this

drift = round([0 cumsum(average_drift(:,2))' ; 0 cumsum(average_drift(:,1))']);
%round([linspace(0,verticalDrift*noPic,noPic);linspace(0,horizontalDrift*noPic,noPic)]);

all_flow = cell2mat(store_flow);
max_deform_x = max(max(all_flow(:,:,1)));
max_deform_y = max(max(all_flow(:,:,2)));
blackPicture = uint8(zeros(picDim(1)+max(drift(1,:)) , picDim(2)+max(drift(2,:)),3));

%am I shifting the right image or am I using the flow on the image before?
%not sure this is working!

for pic = 1:2 %noPic-1
    pngBefore = imread(pngNames(pic+1,:));
    
    pngAfter = blackPicture;
    deform = store_flow{pic};
    
    for ix = 1:picDim(1)
        for iy = 1:picDim(2)  
            %test1 = ix-round(deform(ix,iy,1))+round(max_deform_x)
            %test2 = iy-round(deform(ix,iy,2))+round(max_deform_y)
            pngAfter(ix-round(deform(ix,iy,1))+round(max_deform_x),iy-round(deform(ix,iy,2))+round(max_deform_y),:) = pngBefore(ix,iy,:);                 
        end
    end
    
    imwrite(pngAfter, ['new_' pngNames(pic+1,:)])
end