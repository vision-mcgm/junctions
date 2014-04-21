function [ io,jo,direction ] = getFullFromNeighbourhood2( icord,jcord,direction )
%The function considers each of the pixels surrounding the current position to
% determine where the pixel forming part of the junction is. It first
% consider the pixel along the current direction of travel, then the pixels
% on either side.


global g;

%Parameters

io=0;
jo=0;

found=0;


ringIs=[-1 -1 0 1 1 1 0 -1];
ringJs=[0 1 1 1 0 -1 -1 -1];


ringOffsets=[0 -1 1 -2 2 -3 3 -4];
direction_ring = [5 6 7 8 1 2 3 4 5 6 7 8 1 2 3];

nOffsets=size(ringOffsets,2);

for offsetIndex=1:nOffsets
    if ~found
        thisDirection=direction_ring(direction+ringOffsets(offsetIndex)+4);
        thisI = icord+ringIs(thisDirection);
        thisJ = jcord+ringJs(thisDirection);
        %visited = g.visited(thisI,thisJ,:)
        
        if thisI >0 && thisJ>0 && thisI<g.sizeIm(1) && thisJ < g.sizeIm(2)
            if g.allFull(thisI,thisJ)==1  && ~g.visited(thisI,thisJ,:) && ~g.vertices(thisI,thisJ)
            found=1;
            io=thisI;
            jo=thisJ;
            direction = thisDirection;
            end
        end
    end
end
