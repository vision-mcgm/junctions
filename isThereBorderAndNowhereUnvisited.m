function [ borderTerminator ] = isThereBorderAndNowhereToGo( i,j )
%The function considers the 8 pixels surrounding the current position and
%evaluates whether, for the coordinates of those pixels, there is a border
%of the image. If so, the current position is moved one pixel away from the
%border along the junction and the direction is set to the opposite
%direction of the border

global g;
global q;



border = 0;
somewhereToGo=0;
currI = i;
currJ = j;


for ii=-1:1
    
    for ij=-1:1
        
        if i+ii == 1
            %g.visited(i+ii,j,:) = 1;
            border = 1;
            %currI = 2;
            direction = 5;
            
        elseif i+ii == g.sizeIm(1)
            %g.visited(i+ii,j,:) = 1;
            border = 1;
            %currI = sizeIm(1)-1;
            direction = 1;
        elseif j+ij == 1
            %g.visited(i,j+ij,:) = 1;
            border = 1;
            %currJ = 2;
            direction = 3;
        elseif j+ij == g.sizeIm(2)
            %g.visited(i,j+ij,:) = 1;
            border = 1;
            %currJ = sizeIm(2)-1;
            direction = 7;
        end
        
        if i+ii >0 && j+ij>0 && i+ii <g.sizeIm(1) && j+ij < g.sizeIm(2)
            if ~g.visited(i+ii,j+ij)
                somewhereToGo=1;
            end
        end
        
    end
    
end

borderTerminator=border && ~somewhereToGo;



