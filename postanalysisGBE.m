close all
clear all
clc 

%Load object array of junctions obtained by running findEdges.m 
%load('junctions.mat')
%load('junctions_wildtype.mat')
%load('junctions_GBE_240713_z4.mat')
load('junctions_GBE_OFA_post_1.mat')
load('flow_code/OFA_GBE.mat')

%%%%%%%%%%%%%%%Find the midpoint and length of each junction%%%%%%%%%%%%%%%

%Preallocate array to store the length of junctions
length_junctions = zeros(numberImage,1000);

for time = 1:numberImage
    
    display(time)
    
    %For each junction, generate an image of the single junction. Since the
    %junctions are skeletonized, 'regionprops' can be used to find the
    %perimeter of each junction segment. This is then divided by 2 to get
    %the length
    %(http://blogs.mathworks.com/pick/2012/04/27/calculating-arclengths-made-easy/)
    %Note that the length should always be larger than the pixel length
    %(number of pixels). Also note that the list of junction coordinates,
    %is not necessarily in order.
    
    for n = 1:numberJunctions(time)
        im_length = zeros(sizeImage(time,1),sizeImage(time,2));
        perim = 0;
        length_junction = 0;
        
        %Generate image of single junction
        for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
            im_length(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2)) = 1;
        end
        
        %Only calculate lengths of junctions with at least 1 datapoint
        numCor = size(junctions_array(time,n).junctionCoordinates);
        if numCor(1) > 2
            
            %Calculate length of junction           
            perim = regionprops(im_length,'Perimeter');
            length_junction = perim.Perimeter/2;           
            length_junctions(time,n) = length_junction;
            

            %Calculate midpoint of junction by finding the center of area of pixels
            %in the image and then determining the nearest point that is a
            %pixel. This method will be robust as long as the junction is
            %not excessively curved (this could be made more rigourous).           
            stat = regionprops(im_length,'centroid');            
            %Swap x- and y- values
            swapxy = [stat.Centroid(2),stat.Centroid(1)];
            %Find the point in the list of junction coordinates that is
            %closest to the center of area of the image
            kindex = dsearchn(junctions_array(time,n).junctionCoordinates,swapxy);
            midpointJunction = junctions_array(time,n).junctionCoordinates(kindex,:);

            %For imaging midpoints, uncomment the following
            %imshow(im_length); hold on;
            %for x = 1: numel(stat)
            %plot(midpointJunction(2),midpointJunction(1),'ro');
            %end
            
            %Write length and midpoint into object properties for junction
            junctions_array(time,n).junctionLength = length_junction;
            junctions_array(time,n).midpoint = midpointJunction;
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Assign Junction IDs across time frames%%%%%%%%%%%%%%%%%%%%%%

%Track the midpoint of junctions across time frames and use this to assign
%a unique ID to each junction


%Determine number of junctions in total across frames, as this is the maximum
%number of IDs needed.
lengthlist = 0;
for time = 1:numberImage
    lengthlist = lengthlist + numberJunctions(time);
end

%For the tracking, use the algorithm from http://physics.georgetown.edu/matlab/code.html

%Rewrite data in format that the 'track' code can use as input, positionlistJ(1) and (2) are the coordinates 
%of the midpoint of a junction and (3) is the time.
positionlistJ = zeros(lengthlist,3);

startpoint = 1;
for time = 1:numberImage
    for n = 1:numberJunctions(time)
        if ~isempty(junctions_array(time,n).midpoint),
            positionlistJ(startpoint,1) = junctions_array(time,n).midpoint(1);
            positionlistJ(startpoint,2) = junctions_array(time,n).midpoint(2);
            positionlistJ(startpoint,3) = time;
            startpoint = startpoint + 1;
        end
    end
end

%remove zero entries (should this be done?!) Should it be done with column
%3 instead
positionlistJ(all(positionlistJ==0,2),:)=[];


%Extra warping code

for i=2:numberImage
    uv = estimate_flow_interface(im1,im2,'classic+nl-fast');
    [x,y]=meshgrid(1:sj,1:si);
    points1=...
    for p=1:nPoints
    %Assume Matlab coords
    px=points1(p,1);
    py=points1(p,2);
    newPoints(p,1)=points1(p,1)+uv(px,py,1);
    newPoints(p,2)=points1(p,2)+uv(px,py,2);
    end
    trackerList=... %Add times
    junctionpositions = track(trackerList,8,param);
end

    
    

%Warp midpoints according to flow/deformation matrix found from the OFA

% 
% deformation = 0;
% deformation_x = 0;
% deformation_y = 0;
% 
% for i = 1:numberImage - 1
%     %make deformation matrix cumulative - i.e. image 2 is warped by flow field between image 1 and 2, but image 3 is warped by flow field between image 2 and 3 as well as image 1 and 2, etc. Such that all images are warped 'back' to image 1. 
%     deformation = store_flow{i}-average_drift(i);
%     deformation_x = deformation_x + deformation(:,:,1);
%     deformation_y = deformation_y + deformation(:,:,2);
%     
%     %Note, extracting scattered elements in a matrix requires computing
%     %linear indices of the Matrix. We want the scattered elements
%     %(positionlistJ(find(positionlistJ(:,3) ==
%     %i),1),positionlistJ(find(positionlistJ(:,3) == i),2)) in the matrix
%     %deformation_x. To do this we use the function sub2ind to convert from
%     %row and column subcripts to linear indices.
%     idx = sub2ind(size(deformation_x),positionlistJ(find(positionlistJ(:,3) == i),1) , positionlistJ(find(positionlistJ(:,3) == i),2));
%     positionlistJ(find(positionlistJ(:,3) == i),1) = positionlistJ(find(positionlistJ(:,3) == i),1) - deformation_x(idx);   
%     positionlistJ(find(positionlistJ(:,3) == i),2) = positionlistJ(find(positionlistJ(:,3) == i),2) - deformation_y(idx);
%     
% end
% 
% %Shift all coordinates to positive values
% shift_x = abs(min(positionlistJ(:,1)));
% shift_y = abs(min(positionlistJ(:,2)));
% positionlistJ(:,1) = positionlistJ(:,1) + abs(min(positionlistJ(:,1)));
% positionlistJ(:,2) = positionlistJ(:,2) + abs(min(positionlistJ(:,2)));
% 
% 
% xlim([0 300])
% ylim([0 600])
% for i = 1:82
% pause(0.5)
% scatter(positionlistJ(find(positionlistJ(:,3) == i),1),positionlistJ(find(positionlistJ(:,3) == i),2))
% end

%Track position of midpoints. The second paramater is important as it is 
%the distance beyond which it is assumed that a point in the next time frame 
%is not equal to the point being tracked. When set too small, the same 
%midpoint is given a new track between frames. When set too large, the 
%combinatorics are computationally problematic.

%the memory parameter ensures that a junction undergoing a T1 is assigned
%the same ID before and after the transition, even if there is no junction
%for up to 24 image frames

dd = length(positionlistJ(1,:));
%for mem, generally use 24, but less for Natalia's files
param = struct('mem',6,'dim',dd-1,'good',0,'quiet',0);

junctionpositions = track(positionlistJ,8,param);

%Output is: junctionpositions(1,2) are the midpoint coordinates, (3) is
%time and (4) is a unique ID


%Unwarp midpoints according to flow/deformation matrix found from the OFA

deformation = 0;
deformation_x = 0;
deformation_y = 0;

for i = 1:numberImage - 1
    %make deformation matrix cumulative - i.e. image 2 is warped by flow field between image 1 and 2, but image 3 is warped by flow field between image 2 and 3 as well as image 1 and 2, etc. Such that all images are warped 'back' to image 1. 
    deformation = store_flow{i}-average_drift(i);
    deformation_x = deformation_x + deformation(:,:,1);
    deformation_y = deformation_y + deformation(:,:,2);
    
    %Note, extracting scattered elements in a matrix requires computing
    %linear indices of the Matrix. We want the scattered elements
    %(positionlistJ(find(positionlistJ(:,3) ==
    %i),1),positionlistJ(find(positionlistJ(:,3) == i),2)) in the matrix
    %deformation_x. To do this we use the function sub2ind to convert from
    %row and column subcripts to linear indices.
    
    idx = sub2ind(size(deformation_x),positionlistJ(find(positionlistJ(:,3) == i),1) , positionlistJ(find(positionlistJ(:,3) == i),2));
    positionlistJ(find(positionlistJ(:,3) == i),1) = positionlistJ(find(positionlistJ(:,3) == i),1) + deformation_x(idx);   
    positionlistJ(find(positionlistJ(:,3) == i),2) = positionlistJ(find(positionlistJ(:,3) == i),2) + deformation_y(idx);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Test how many frames each midpoint is tracked for (this tells how long
%each junction persists - most should persist for the duration of the
%movie)
for k = 1:max(unique(junctionpositions(:,4)))
    count(k) = length(find(junctionpositions(:,4)==k));
end


%%%%%%%%%%%%%%%Assign IDs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Assign each junction an individual ID based on tracking the midpoints. In
%the array of junctions, each junction is identified by the time and an
%index n. This index can change between time frames and does therefore not
%track the junctions over time.

for time = 1:numberImage
    
    %Select rows from junctionpositions that correspond to this time
    sub_junctionpositions = junctionpositions(junctionpositions(:,3) == time,:);
    
    for n = 1:length(junctions_array)  %numberJunctions*2
        %Match midpoint coordinates for junction with entry in
        %junctionpositions
        if ~isempty(junctions_array(time,n).midpoint),
            index_track = mfind(sub_junctionpositions(:,1:2),junctions_array(time,n).midpoint);            
            %Assign ID to junction objects.
            junctions_array(time,n).junctionID = sub_junctionpositions(index_track,4);
            %Create list that allows translation from n index to ID index 
            %and back for each time point
            id_of_n(time,n) = sub_junctionpositions(index_track,4);
            n_of_id(time,sub_junctionpositions(index_track,4)) = n;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%Generate list of junction lengths according to IDs%%%%%%%%%%

list_ids = unique(id_of_n)';
list_ids(list_ids == 0) = [];
length_array = zeros(numberImage,length(list_ids));
        
for ids = list_ids
            for t = 1:numberImage
                if ~(n_of_id(t,ids) == 0)
                length_array(t,ids) = length_junctions(t,n_of_id(t,ids));               
                end
            end         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Find neighbours of junctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preallocate array of neighbours. The value will be 1 if n,m are first neighbours
%and 0 if they're not.
junction_neighbours = zeros(numberImage,max(numberJunctions),max(numberJunctions));

for time = 1:numberImage;
    display(time)

   for m = 1:numberJunctions(time)          
       for n = 1:numberJunctions(time);
            
       neighbour = 0;
       
       %Do not count junctions as having themselves as neighbours
       if n == m         
        neighbour = 0;
       else                     
       if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,m).vertex1)  
             if junctions_array(time,n).vertex1 == junctions_array(time,m).vertex1
             neighbour = 1;
             end
       end               
       if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,m).vertex2)
            if junctions_array(time,n).vertex1 == junctions_array(time,m).vertex2
            neighbour = 1;
            end
       end
       if ~isempty(junctions_array(time,n).vertex2) && ~isempty(junctions_array(time,m).vertex1)
            if junctions_array(time,n).vertex2 == junctions_array(time,m).vertex1
            neighbour = 1;
            end
       end
       if ~isempty(junctions_array(time,n).vertex2) && ~isempty(junctions_array(time,m).vertex2)
            if junctions_array(time,n).vertex2 == junctions_array(time,m).vertex2
            neighbour = 1;
            end
       end        
       end
       
       if neighbour == 1,
         junction_neighbours(time,m,n) = 1;  
       end          
       
       end         
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('junctions_GBE_OFA_post_1.mat')
%save('junctions_wildtype_post.mat')
%save('junctions_GBE_240713_z4_data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%