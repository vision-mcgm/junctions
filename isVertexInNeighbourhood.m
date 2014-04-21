function [ tf,iver,jver ] = isVertexInNeighbourhood( i,j )
%The function considers the 8 pixels surrounding the current position and
%evaluates whether, for the coordinates of those pixels, the entries in
%g.vertices have the value 1. If so, it returns tf = 1, indicating that
%there is a vertex in the neighbourhood, as well as the coordinates of the
%vertices.

global g;

tf=0;
iver = 0;
jver = 0;

for ii=-1:1
    for ij=-1:1
          if i+ii >0 && j+ij>0 && i+ii <g.sizeIm(1) && j+ij < g.sizeIm(2)
            if g.vertices(i+ii,j+ij)==1
                
                tf=1;
                iver = i+ii;
                jver = j+ij;
                
            end
        end
    end
    
end

