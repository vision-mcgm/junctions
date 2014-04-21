function [ junctions , vertexcoordinates , numJunctions ] = findEdges(im) % insert im back for calling from image_stack_analysis
%Find edges

%vertexcoordinates=0;
global g;
global q;
q=5;


%Load first image of movie
%im = imread('181012 SPE2 bd_fate files/bd_fate001.png');  %imread('bd_fate00-04.tif',1); %comment out when using image_stack_analysis

%Define subset of image to test code for
%sub=im(100:400,100:400,:);
%im = sub;

%Parameters
workIm=im;
sizeIm=size(workIm);
g.sizeIm=sizeIm;


[si sj c]=size(workIm);
g.allFull=zeros(si,sj);

%Define and display the green colour channel, contains both vertices and a
%few of the junctions coloured green/yellow in the original image
gch=workIm(:,:,2);
rch=workIm(:,:,1);
bch=workIm(:,:,3);

%The vertices are the entries in the green colour channel given by
%intensity 255, define coordinates for the vertices
%[iVg jVg]=find(gch==255); %finds indices for entries in green colour channel of intensity 255

%Find which pixels are in all
verticesi=[];
verticesj=[];
g.vertices=zeros(si,sj);
for ii=1:si
    for ij=1:sj
        if workIm(ii,ij,1)==255 && workIm(ii,ij,2)==255 && workIm(ii,ij,3)==255
            verticesi=[verticesi ii];
            verticesj=[verticesj ij];
            g.vertices(ii,ij,:)=1;
        end
    end
end

vertexid = zeros(length(verticesi),1);
vertexcoordinates = [verticesi' verticesj' vertexid];


%Code


%[iJ jJ]=find(rch==255);

iJ=[];
jJ=[];
%tic;
for ii=1:si
    for ij=1:sj
        if ~(im(ii,ij,1)==0 && im(ii,ij,2)==0 && im(ii,ij,3)==0)
            iJ=[iJ ii];
            jJ=[jJ ij];
        end
    end
end

%toc
iJ = iJ';
jJ = jJ';
sizeiJ=size(iJ,1);

for n=1:sizeiJ
g.allFull(iJ(n),jJ(n),:)=1;
end


%Only choose starting positions from list of coordinates not containing
%vertex coordinates. Create iJnv and jJnv to be the lists without vertex
%coordinates

junctioncoords = [iJ jJ];
vertexcoords = [verticesi' verticesj'];

[~, ia, ~] = intersect(junctioncoords,vertexcoords,'rows');
junctioncoords(ia,:) = [];
iJnv = junctioncoords(:,1);
jJnv = junctioncoords(:,2);


%n=size(iJ,1);

%[si sj c]=size(im);

empty=0;

g.visited=zeros(si,sj);


%Initialize array of objects - For objects the pre-allocation works by assigning one of the objects to the very last field in the array. Matlab then fills the other fields before that with objects (handles) that it creates by calling the constructor of that object with no arguments.
junctions(1,1000) = junction();
for n = 1:1000
junctions(n) = junction();
end

%Loop over junctions
    
ctjunction = 0;
numJunctions = 0;
alldone = 0;

for ctjunction = 1:1000
   
    
    index = 0;
    for index = 1:length(iJnv)         
         startI=iJnv(index);
         startJ=jJnv(index); 
         direction =3;
         if g.visited(startI,startJ)==0,break,
         elseif index == length(iJnv), 
             numJunctions = ctjunction-1;
             alldone = 1;
         end
    end
    

    %If the last iteration found that no points were left unvisited, break
    if alldone == 1, break, end
     
currI = startI;
currJ = startJ;


finished=0;
border = 0;
vertexfound=0;
it = 1;
terminators=0;


while ~finished
    %pause(0.005) %This is just to give Matlab time to update the image with imshow 
    
    
    %store junction coordinates
    junctions(ctjunction).junctionCoordinates(it,1) = currI;
    junctions(ctjunction).junctionCoordinates(it,2) = currJ;
    
    %add junction coordinate to visited
    g.visited(currI,currJ) = 1;
    workIm(currI,currJ,:)=[0 255 0];

    %plots detected pixels as we go along
%     if mod(it,1)==0
%     imshow(imCross(workIm,currJ,currI),'InitialMagnification',100);
%     end
        
    if isThereBorderAndNowhereUnvisited(currI,currJ)
        terminators=terminators+1;
    end
    
    %oldBorder=border;
    %[border,currI_dontcare,currJ_dontcare,direction] = isThereBorder(currI,currJ,direction,sizeIm,border);
    %if border >oldBorder
     %   borderFound=1;
    %else borderFound=0;
    %end
    

    [tf,iver,jver] = isVertexInNeighbourhood(currI,currJ);
    if tf == 1,
        if vertexfound == 0,
           vertexfound = vertexfound+1;           
           junctions(ctjunction).vertex1(1) = iver;
           junctions(ctjunction).vertex1(2) = jver;
           
           %check whether vertex has been added to list of junction
           %coordinates, if not, add it
           
           if ~(junctions(ctjunction).junctionCoordinates(it,1) == iver && junctions(ctjunction).junctionCoordinates(it,2) == jver),
               it = it+1;
               junctions(ctjunction).junctionCoordinates(it,1) = iver;
               junctions(ctjunction).junctionCoordinates(it,2) = jver;
           end
           
           %if it's the first vertex, jump back to starting point and walk
           %in the opposite direction along junction.
           if it>1 %or should that be one? 
           currI = startI;
           currJ = startJ;
           direction = 7;
           %change the order of pixels such that they are in order
           %rather than with a jump
           junctions(ctjunction).junctionCoordinates = junctions(ctjunction).junctionCoordinates(end:-1:1,:);
           end
        else
           %check that vertex is not the one already found
           if ~(iver==junctions(ctjunction).vertex1(1) && jver==junctions(ctjunction).vertex1(2))        
           vertexfound = 2;
           junctions(ctjunction).vertex2(1) = iver;
           junctions(ctjunction).vertex2(2) = jver;
           
           %check whether vertex has been added to list of junction
           %coordinates, if not, add it
           
           if ~(junctions(ctjunction).junctionCoordinates(it,1) == iver && junctions(ctjunction).junctionCoordinates(it,2) == jver),
               it = it+1;
               junctions(ctjunction).junctionCoordinates(it,1) = iver;
               junctions(ctjunction).junctionCoordinates(it,2) = jver;
           end
           
           end
        end
    end
    
    if vertexfound == 2 || vertexfound+terminators == 2 || terminators==2
        finished = 1;
        break,
    end
     
    [nextI, nextJ, direction] = getFullFromNeighbourhood2(currI,currJ,direction);
    
    if nextI==0 || nextJ ==0
        break %Can't find a next pixel. Break the loop. 
    end
    
    currI=nextI;
    currJ=nextJ;
    
    it = it+1;
     
end

    
end
end
        
        
        