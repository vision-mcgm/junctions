function [ im ] = imCross( im,x,y )
%imCross(im,x,y)
x=round(x);
y=round(y);

d=1;

is=size(im,1);
js=size(im,2);


if x > js || x <2 || y >is || y <2
    
im(y:y,x,:)=255;
im(y,x:x,:)=255;
else
im(y-d:y+d,x,:)=255;
im(y,x-d:x+d,:)=255;
end

% im(y-d:y+d,x,1)=min(255,im(y-d:y+d,x,1)+10);
% im(y,x-d:x+d,1)=min(255,im(y,x-d:x+d,1)+10);

end

