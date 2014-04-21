%im1 = imread('image1.png');
%im2 = imread('image2.png');
%uv = estimate_flow_interface(im1,im2,'classic+nl-fast');

% Display estimated flow fields
%figure; subplot(1,2,1);imshow(uint8(flowToColor(uv))); title('Middlebury color coding');
%subplot(1,2,2); plotflow(uv);   title('Vector plot');

close all
clear all

store_flow = cell(81,1);
average_drift = zeros(81,2);

for i = 1:81
    fname = 'bd_fate%.3d.png';
    clearvars im1 im2 uv
    im1 = imread(sprintf(fname,i));
    im2 = imread(sprintf(fname,i+1));
    uv = estimate_flow_interface(im1,im2,'classic+nl-fast');
    store_flow{i} = uv;
    average_drift(i,:) = mean(mean(uv));   
end

for time = 1:81
    
    filename = ['flow', num2str(time)];
    clear flow
    flow = store_flow{time}-average_drift(time);
    pic = figure;
    plotflow(flow);
    print(pic,'-dpng',filename);
    
end


