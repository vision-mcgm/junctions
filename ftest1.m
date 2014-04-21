function [  ] = ftest1( im1,im2,store_flow)
%Test flow field in various ways

%uv = estimate_flow_interface(im1,im2,'classic+nl-fast');

%Reverse warp

testfield = store_flow;
im1=im2double(im1);
im2=im2double(im2);
im1=rgb2gray(im1);
im2=rgb2gray(im2);

[si sj]=size(im1);

% xf=testfield(:,:,1);
% yf=testfield(:,:,2);

xf=ones(si,sj)*100;
yf=ones(si,sj)*100;

[x,y]=meshgrid(1:sj,1:si);

[qx,qy]=meshgrid(1:sj,1:si);

%R1=interp2(im2,x-testfield(:,:,1),y-testfield(:,:,2));
%R1=interp2(im2,x-xf,y-yf); %Apply xf and yf to im2 and interpolate.
W1=interp2(x,y,im1,qx,qy);
figure(1)
imshow(im1)
figure(2)
imshow(W1)
keyboard




end

