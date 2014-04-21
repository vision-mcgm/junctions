close all
clear all
clc
tic

%Load object array of junctions obtained by running findEdges.m
%load('junctions.mat')
%load('junctions_wildtype.mat')
%load('junctions_GBE_240713_z4.mat')
%load('flow_code/OFA_GBE_new.mat')
ffolder = '181012 SPE2 bd_fate files';
loadname = [ffolder '_junctions' '.mat'];
load(loadname);

%%%%%%%%%%%%%%%Find the midpoint and length and orientation of each junction%%%%%%%%%%%%%%%

for time = 1:numberImage
    
    display(time)
    
    for j = 1:numberJunctions(time)
              
        %Only calculate lengths of junctions with at least 1 datapoint
        if size(junctions_array(time,j).junctionCoordinates,1) > 2
            
            %Calculate length of junction
            length_junction = 0;
            coordinates = junctions_array(time,j).junctionCoordinates;
            pixelInSamePoleAsBefore = mod(diff(sum(coordinates,2)),2); %1 if this pixel is in the same pole. 0 if new pole started
            
            poleLength=1; %a pole is a block of pixels that are connected (not diagonal)
            for i=1:length(pixelInSamePoleAsBefore)
                if(pixelInSamePoleAsBefore(i)) 
                    poleLength = poleLength+1; %make this pole 1 pixel longer
                else
                    length_junction = length_junction + sqrt(poleLength^2+1); %add diagonal length of pole to junction length
                    poleLength = 1;
                end
            end
            length_junction = length_junction + sqrt(poleLength^2+1);
            
            %save junction length
            junctions_array(time,j).junctionLength = length_junction;
            
            
            %Calculate midpoint of junction as central pixel
            junctions_array(time,j).midpoint = coordinates(ceil(end/2),:);
            
            
            %Calculate orientation of junction if both vertices exist
            if ~isempty(junctions_array(time,j).vertex1) && ~isempty(junctions_array(time,j).vertex2)
                deltax = junctions_array(time,j).vertex2(1)-junctions_array(time,j).vertex1(1);
                deltay = junctions_array(time,j).vertex2(2)-junctions_array(time,j).vertex1(2);

                if(deltax==0)
                    junctions_array(time,j).angle = pi;
                else
                    junctions_array(time,j).angle = atan2(deltay,-deltax);
                end
            end
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Assign Junction IDs across time frames%%%%%%%%%%%%%%%%%%%%%%

%Track the midpoint of junctions across time frames and use this to assign
%a unique ID to each junction


%Determine number of junctions in total across frames, as this is the maximum
%number of IDs needed.
lengthlist = sum(numberJunctions);

%For the tracking, use the algorithm from http://physics.georgetown.edu/matlab/code.html

%Rewrite data in format that the 'track' code can use as input, positionlistJ(1) and (2) are the coordinates
%of the midpoint of a junction and (3) is the time.
positionlistJ = zeros(lengthlist,3);

startpoint = 1;
for time = 1:numberImage
    for j = 1:numberJunctions(time)
        if ~isempty(junctions_array(time,j).midpoint),
            positionlistJ(startpoint,1) = junctions_array(time,j).midpoint(1);
            positionlistJ(startpoint,2) = junctions_array(time,j).midpoint(2);
            positionlistJ(startpoint,3) = time;
            startpoint = startpoint + 1;
        end
    end
end

%remove zero entries (should this be done?!) Should it be done with column
%3 instead
positionlistJ(all(positionlistJ==0,2),:)=[];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Test how many frames each midpoint is tracked for (this tells how long
%each junction persists - most should persist for the duration of the
%movie)
% for k = 1:max(unique(junctionpositions(:,4)))
%     count(k) = length(find(junctionpositions(:,4)==k));
% end


%%%%%%%%%%%%%%%Assign IDs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Assign each junction an individual ID based on tracking the midpoints. In
%the array of junctions, each junction is identified by the time and an
%index n. This index can change between time frames and does therefore not
%track the junctions over time.
id_of_n = zeros(numberImage,max(numberJunctions));
n_of_id = zeros(numberImage,max(junctionpositions(:,4)));

for time = 1:numberImage
    
    %Select rows from junctionpositions that correspond to this time
    sub_junctionpositions = junctionpositions(junctionpositions(:,3) == time,:);
    
    for j = 1:numberJunctions(time) 
        
        %Match midpoint coordinates for junction with entry in junctionpositions
        if ~isempty(junctions_array(time,j).midpoint),
            index_track = mfind(sub_junctionpositions(:,1:2),junctions_array(time,j).midpoint);
            %Assign ID to junction objects.
            junctions_array(time,j).junctionID = sub_junctionpositions(index_track,4);
            %Create list that allows translation from n index to ID index
            %and back for each time point
            id_of_n(time,j) = sub_junctionpositions(index_track,4);
            n_of_id(time,sub_junctionpositions(index_track,4)) = j;
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
            length_array(t,ids) = junctions_array(t,n_of_id(t,ids)).junctionLength;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Find neighbours of junctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Preallocate cell array of neighbours. junction_neighbours(time,junction)
%has an arry with the neighbours of that junction at that time.

for time = 1:numberImage;
    fprintf(['finding neighbours. time = ' num2str(time) '\n'])
    
    for j1 = 1:numberJunctions(time) 
        for j2 = j1+1:numberJunctions(time); %only loop over different j1 and j2
            
            neighbour = false;
            
            if ~isempty(junctions_array(time,j2).vertex1) 
                if ~isempty(junctions_array(time,j1).vertex1) && all(junctions_array(time,j2).vertex1 == junctions_array(time,j1).vertex1)
                    neighbour = true;
                end
                if ~isempty(junctions_array(time,j1).vertex2) && all(junctions_array(time,j2).vertex1 == junctions_array(time,j1).vertex2)
                	neighbour = true;
                end
            end
            
            if ~isempty(junctions_array(time,j2).vertex2) 
                if ~isempty(junctions_array(time,j1).vertex1) && all(junctions_array(time,j2).vertex2 == junctions_array(time,j1).vertex1)
                    neighbour = true;
                end
                
                if ~isempty(junctions_array(time,j1).vertex2) && all(junctions_array(time,j2).vertex2 == junctions_array(time,j1).vertex2)
                    neighbour = true;
                end
            end
               
            if neighbour
                junctions_array(time,j1).neighboursID = [junctions_array(time,j1).neighboursID junctions_array(time,j2).junctionID];
                junctions_array(time,j2).neighboursID = [junctions_array(time,j2).neighboursID junctions_array(time,j1).junctionID];
            end
            
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

savename = [ffolder '_junctionsWithProps' '.mat'];
save(savename,'junctions_array','vertex_array','numberJunctions','numberImage','id_of_n','n_of_id','length_array')

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
