clear all
close all

load('junction_data.mat')


%%%%correlation analysis

%Choosing junction to compare to
%user selects approximate midpoint of junction. Closest midpoint is then
%found and the junction is chosen

%Plot junction at time = 1 to choose from
time = 1;
image_junctions = zeros(sizeImage(time,1),sizeImage(time,2));

for n = 1:numberJunctions(time)
        %assign colour
        C = 100;
        
        %Colour image of single junction
        for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
            image_junctions(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2),:) = C;
        end
        
end
           
cmp = colormap(jet(100)); %hsv
cmp(1,:) = [0 0 0];
imshow(image_junctions,cmp);


[xm ym] = ginput(1)
chosen_point = [xm,ym];
%hold on
%plot(xm,ym,'ms')
chosen_point_swap = [ym,xm];

%find junction midpoint closest to the point chosen
sub_junctionpositions = junctionpositions(junctionpositions(:,3) == time,:);
index_p = dsearchn(sub_junctionpositions(:,1:2),chosen_point_swap);
chosen_junction_id = sub_junctionpositions(index_p,4);
%hold on
%plot(junctions_array(time,n_of_id(time,chosen_junction_id)).midpoint(2),junctions_array(time,n_of_id(time,chosen_junction_id)).midpoint(1),'b*')

%%%%%%%%%%%%%%%%%%%%%
%Plot chosen junction in white. Calculate correlation coefficients for all
%other junctions.

chosen_junction_id = 1588 %358(184);

window_size = 25 %input('please enter window_size \n');     
chosen_data = length_array(:,chosen_junction_id);

%include colorbar in separate window
%cmp = colormap(jet(103)); %hsv
%ch = colorbar('YLim',[2 102],...                        &# The axis limits
%              'YTick',[2 51 102],...                    %# The tick locations
%              'YTickLabel',{'-1, anti-correlation','0, no correlation','+1, correlation'});  %# The tick labels
 
clear C_vector
%figure(1)
%imshow('colorbar.png')
%figure(2)
for time = 1:numberImage-window_size
    %figure()
    image_junctions = zeros(sizeImage(time,1),sizeImage(time,2));
    
    %The option 'pic' and 'saveas' are if one wants to create a movie out of
    %the images
    filename = ['image', num2str(time)];
    %pic=figure();
    
 
    
    for n = 1:numberJunctions(time)

        %Determine correlation between this junction and the chosen
        %junction
        
%         for n = 1:numberJunctions(time)
%         for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
%             image_junctions(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2),:) = 102;
%         end   
%         end
        
        
        %Colour image of single junction
        %for n = n_of_id(time,[17 21 24 28])

        if ~(id_of_n(time,n) == 0)
        other_data = length_array(:,id_of_n(time,n)); 
        
        %window the data
        chosen_data_w = chosen_data(time:time+window_size);
        other_data_w = other_data(time:time+window_size);
        
        [C,lags] = xcorr(chosen_data_w-mean(chosen_data_w),other_data_w-mean(other_data_w),0,'coeff');
        %ind = [2:1:101];
        C_vector(time,n) = C;
        C = 100*(C+1)/2+2;
        %else
            %C = 102;
            
        for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
            image_junctions(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2),:) = C;
        end
        end
        
        %end
        %should the time be at the midpoint of the window, at the start, or
        %at the end?

        % for i = 1:length(raw_data(:,2))-window_size
        % [c,lags] = xcorr(chosen_data(i:i+window_size)-mean(chosen_data(i:i+window_size)),other_data(i:i+window_size)-mean(other_data(i:i+window_size)),0,'coeff');
        % c_vector(i) = c;
        % end


        if ~(n_of_id(time,chosen_junction_id) == 0)
        %Colour chosen junction
        for i = 1:length(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(:,1))
            image_junctions(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,1),junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,2),:) = 103;
        end
        end
    end
        
        
    %end
    
    %maxC = max(max(image_junctions));
    
    %caxis([-101,101])
    cmp = colormap(jet(103)); %hsv
    cmp(1,:) = [0 0 0];
    cmp(103,:) = [1 1 1];
    cmp(102,:) = [0.5 0.5 0.5];
    image_junctions_sub = image_junctions(170:250,240:320,:);
    pic = imshow(image_junctions_sub,cmp,'InitialMagnification',200); %colorbar;
    saveas(pic,filename,'png');
    
end
