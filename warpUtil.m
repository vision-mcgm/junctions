function [ y ] = warpUtil( im,vx,vy )
C= imread('pout.tif'); % test image
[x, y] = meshgrid(1:size(im,2), 1:size(im,1));



% compute the warped image - the subtractions are because we're specifying
% where in the original image each pixel in the new image comes from
y = interp2(double(im), x-vx, y-vy);

% display the result



end

