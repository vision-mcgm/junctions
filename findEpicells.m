clear all
close all
clc
tic

%load junctions
ffolder = '181012 SPE2 bd_fate files'; %'GBE_OFA_stable';%'bdfate GBE 240713 z4';%'181012 SPE2 bd_fate files';%'bd_fate_high_res'; %'DECADGFP ROKiR 240113 SPE2 bd_fate'; %'bd_fate_T1' to use other folder
loadname = [ffolder '_junctionsWithProps' '.mat'];
load(loadname)

numberEpicellsMany = 1000;
numberEpicells = zeros(numberImage,1);

%Initialize array of objects - For objects the pre-allocation works by assigning one of the objects to the very last field in the array. Matlab then fills the other fields before that with objects (handles) that it creates by calling the constructor of that object with no arguments.
epicells_array(numberImage,numberEpicellsMany) = epicell();
for time = 1:numberImage
    for n = 1:numberEpicellsMany
        epicells_array(time,n) = epicell();
    end
end


for time = 1:numberImage
    cellNumber = 1;

    %for each junction, go clockwise and counterclockwise
    for j = 1:numberJunctions(time)
        if ~isempty(junctions_array(time,j).angle) %only start if both vertices exist - ie if angle exists
            
            for rotation = 1:2
                junctionIDInCell = id_of_n(time,j);
                
                %the two last vertices in 'verticesInCell determines directions you walk - next junction is always the neighbour linked to the last vertex that is most to the right when coming from the second last vertex
                if(rotation ==1)
                    verticesInCell = junctions_array(time,j).vertex1;
                    verticesInCell(2,:) = junctions_array(time,j).vertex2;
                    presentAngle = junctions_array(time,j).angle;
                else
                    verticesInCell = junctions_array(time,j).vertex2;
                    verticesInCell(2,:) = junctions_array(time,j).vertex1;
                    presentAngle = junctions_array(time,j).angle+pi; %angle is pi different, if you move from vertex 2 to vertex 1
                end
                
                for step = 1:8
                    angleToNextJunction = 10; %the next junction is the one with the smallest angle.
                    nextJunction = [];
                    nextAngle = [];
                    nextVertex = [];
                    
                    
                    for neighbourJunction = junctions_array(time,n_of_id(time,junctionIDInCell(end))).neighboursID %go though neigbours
                        if ~isempty(junctions_array(time,n_of_id(time,neighbourJunction)).angle) %only consider neighbour if both vertices exist - ie if angle exists
                            
                            
                            if(junctions_array(time,n_of_id(time,neighbourJunction)).vertex1 == verticesInCell(end,:)) %if vertex 1 of this neighbour is connected to the present vertex
                                neighbourangle = junctions_array(time,n_of_id(time,neighbourJunction)).angle;
                                angleDiff = mod(3*pi+presentAngle - neighbourangle , 2*pi);
                                
                                %if this is the smallest angle so far
                                if(angleDiff<angleToNextJunction)
                                    angleToNextJunction = angleDiff;
                                    nextJunction = neighbourJunction;
                                    nextAngle = neighbourangle;
                                    nextVertex = junctions_array(time,n_of_id(time,neighbourJunction)).vertex2;
                                end
                                
                            elseif(junctions_array(time,n_of_id(time,neighbourJunction)).vertex2 == verticesInCell(end,:)) %if vertex 2 of this neighbour is connected to the present vertex
                                neighbourangle = junctions_array(time,n_of_id(time,neighbourJunction)).angle+pi; %negative orientation if you go from vertex 2 to 1
                                angleDiff = mod(3*pi+presentAngle - neighbourangle , 2*pi);
                                
                                %if this is the smallest angle so far
                                if(angleDiff<angleToNextJunction)
                                    angleToNextJunction = angleDiff;
                                    nextJunction = neighbourJunction;
                                    nextAngle = neighbourangle;
                                    nextVertex = junctions_array(time,n_of_id(time,neighbourJunction)).vertex1;
                                end
                            end
                        end
                    end
                    
                    if(isempty(nextJunction) || n_of_id(nextJunction)<j) %if next junction does not exist, or if it is less than the starting junction (and therefore already looked at)
                        break;
                    else %if there was a neighbour junction from that vertex
                        junctionIDInCell(end+1) = nextJunction;
                        
                        if(verticesInCell(1,:) == nextVertex) %if the next vertex is the same as the starting vertex, the epicell is found
                            epicells_array(time,cellNumber).junctionIDs = junctionIDInCell;
                            epicells_array(time,cellNumber).vertices = verticesInCell;
                            cellNumber = cellNumber+1;
                            break;
                        else %if we are not back yet
                            verticesInCell(end+1,:) = nextVertex;
                            presentAngle = nextAngle;
                        end
                    end
                end
            end
        end
    end
    numberEpicells(time) = cellNumber -1;
    
end

toc