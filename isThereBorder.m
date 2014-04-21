function [ border,currI,currJ,direction ] = isThereBorder( i,j,direction,sizeIm,border )
%The function considers the 8 pixels surrounding the current position and
%evaluates whether, for the coordinates of those pixels, there is a border
%of the image. If so, the current position is moved one pixel away from the
%border along the junction and the direction is set to the opposite
%direction of the border

global g;

%border = 0;
currI = i;
currJ = j;


for ii=-1:1
    
    for ij=-1:1
        
        if i+ii == 1
            %g.visited(i+ii,j,:) = 1;
            border = 1;
            %currI = 2;
            direction = 5;
        elseif i+ii == sizeIm(1)
            %g.visited(i+ii,j,:) = 1;
            border = 1;
            %currI = sizeIm(1)-1;
            direction = 1;
        elseif j+ij == 1
            %g.visited(i,j+ij,:) = 1;
            border = 1;
            %currJ = 2;
            direction = 3;
        elseif j+ij == sizeIm(2)
            %g.visited(i,j+ij,:) = 1;
            border = 1;
            %currJ = sizeIm(2)-1;
            direction = 7;
        end
      
    end

end

