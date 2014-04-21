%Assign colour according to ID

for time = 1:numberImage
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
        
        C = id_of_n(time,n);
            
        for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
            image_junctions(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2),:) = C;
        end
        
        %end
        %should the time be at the midpoint of the window, at the start, or
        %at the end?

        % for i = 1:length(raw_data(:,2))-window_size
        % [c,lags] = xcorr(chosen_data(i:i+window_size)-mean(chosen_data(i:i+window_size)),other_data(i:i+window_size)-mean(other_data(i:i+window_size)),0,'coeff');
        % c_vector(i) = c;
        % end


    end
        
        
    %end
    
    %maxC = max(max(image_junctions));
    
    %caxis([-101,101])
    cmp = colormap(jet(1600)); %hsv
    cmp(1,:) = [0 0 0];
    %cmp(103,:) = [1 1 1];
    %cmp(102,:) = [0.5 0.5 0.5];
    image_junctions_sub = image_junctions(170:250,240:320,:);
    pic = imshow(image_junctions_sub,cmp,'InitialMagnification',200);
    %pic = imshow(image_junctions,cmp); %colorbar;
    saveas(pic,filename,'png');
    
end