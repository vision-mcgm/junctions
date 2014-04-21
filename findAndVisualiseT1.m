%Code to identify junctions that shrink to a four-way vertex and identify
%junctions that undergo a T1 transition.

%Border junctions are excluded.

%The time at which T1 transitions occur are recorded and the angle of
%junctions as they approach and leave a T1 transition is determined.

%All T1 transitions are then plotted to examine that they indeed are T1
%events.

%The code considers the junctions one at a time (iterates over all
%junctions) and then considers each time point (from time =
%2:numberImage-1). Junctions are considered if they currently ('time') have a
%length>0 and in the next step ('time+1') have lenght=0. The first time
%point at which the junction length increases to length>0 is denotes
%'time_after'.

close all
clear all
%load('junctions_GBE_240713_z4.mat')

%%%%%Identify border junctions (in order to exclude these from subsequent
%%%%%analysis).

border1 = 2;
border2 = 2;
border3 = sizeImage(time,1)-1;
border4 = sizeImage(time,2)-1;

for time = 1:numberImage
    for m = 1:numberJunctions(time)
        if any(junctions_array(time,m).junctionCoordinates(:,1) == border1) || any(junctions_array(time,m).junctionCoordinates(:,1) == border3) || any(junctions_array(time,m).junctionCoordinates(:,2) == border2) || any(junctions_array(time,m).junctionCoordinates(:,2) == border4)
            bordertest(time,m) = 1;
        end
    end
end


%%%%Initialize variables
%countshrink = 0;
%shrink_ids = [];
%shrink_time = [];
count_almost_T1 = 0;
almost_T1_ids = [];
almost_T1_time = [];
countT1 = 0;
T1_ids = [];
T1_time = [];
T1_list = zeros(numberImage,max(numberJunctions));
count_t1 = 1;
total_junctions = [];
angle_before_T1 = [];
sign_before_T1 = [];


%%%%From the data exclude the following junctions: borderjunctions, any
%%%%junctions that don't have two vertices to start with, junctions that
%%%%participate in cell divisions (these jump from being long to being zero
%%%%in one time step).


%Note: are there any non-border junctions that don't have two vertices to start with? Is the requirement that long junctions can't suddenly shrink too
%stringent?

for time = 2:numberImage-1
    
    for m = 1:numberJunctions(time)
        
        %exclude borderjunctions
        if bordertest(time,m) == 0
            
            %only consider junctions that have two vertices to start with
            %(generally redundant - this criteria has not been necessary
            %for any of the movies so far)
            if ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                
                if ~(id_of_n(time,m) == 0)
                    mid = id_of_n(time,m);
                    
                    total_junctions = [total_junctions mid];
                    
                    %exclude junctions involved in cell divisions (by assuming that
                    %a junction can't jump from being long (more than 20) to zero in one time step)
                    if (length_array(time,mid) < 20) && length_array(time-1,mid) < 20
                        
                        %%%%%check whether the junction length shrinks to zero (indicates four-way vertx) in next time interval
                        if length_array(time+1,mid) == 0
                            %countshrink = countshrink + 1;
                            %shrink_ids = [shrink_ids mid];
                            %shrink_time = [shrink_time time+1];
                            
                            
                            %%%%check if the junction then undergoes a T1 transition
                            
                            %find first timepoint (after current time point) where the junction has nonzero length and find the neighbours at this time point
                            nj1 = find(~(length_array(:,mid) == 0));
                            nj2 = nj1(find(nj1 > time));
                            
                            
                            %If there is no such time, the junction has disappeared and
                            %is not a T1 transition.
                            
                            %%%%Determine id of neigbhouring junctions at 'time' and 'time_after'
                            
                            if ~isempty(nj2)
                                time_after = nj2(1);
                                
                                
                                %To exclude cell divisions, require all four neighbours to be the same before and
                                %after. Events that fulfill these criteria (junctions that
                                %contract to a four-way vertex) will be referred to as
                                %'almost_T1'. Record number, id of junction, and time of such events.
                                if isequal(sort(junctions_array(time,m).neighboursID),sort(junctions_array(time_after,n_of_id(time_after,mid)).neighboursID))
                                    count_almost_T1 = count_almost_T1 + 1;
                                    almost_T1_ids = [almost_T1_ids mid];
                                    almost_T1_time = [almost_T1_time time];
                                    
                                    
                                    %For true T1s, the junction itself has the same four
                                    %neighbouring junctions, but these junctions swap
                                    %neighbours. If neighbour 1 and 2 shared a vertex before,
                                    %after the T1 transition 1 will instead share a vertex with
                                    %either neighbour 3 or 4. Check for this by idenfiying the
                                    %neighbour of junction 1 that is also a neighbour of the
                                    %junction undergoing the T1 transition, and checking if
                                    %this is the same after the event.
                                    
                                    neighbour1 = junctions_array(time,m).neighboursID(1);
                                    
                                    %Note, we've already checked that
                                    %junctions_array(time,m).neighboursID = junctions_array(time_after,n_of_id(time_after,mid)).neighboursID)
                                    neighbour_before = 0;
                                    neighbour_after = 0;
                                    neighbour_before_index = 0;
                                    neighbour_after_index = 0;
                                    
                                    neighbours_before_index = find(ismember(junctions_array(time,m).neighboursID,junctions_array(time,n_of_id(time,neighbour1)).neighboursID));
                                    neighbours_after_index = find(ismember(junctions_array(time,m).neighboursID,junctions_array(time_after,n_of_id(time_after,neighbour1)).neighboursID));
                                    neighbour_before = junctions_array(time,n_of_id(time,neighbour1)).neighboursID(neighbours_before_index);
                                    neighbour_after = junctions_array(time_after,n_of_id(time_after,neighbour1)).neighboursID(neighbours_after_index);
                                    
                                    %If there is a neighbour exchange, count it as a T1 and
                                    %record junction ID and time.
                                    if ~(neighbour_before == neighbour_after)
                                        countT1 = countT1 + 1;
                                        T1_ids = [T1_ids mid];
                                        T1_time = [T1_time time];
                                        
                                        
                                        %%%%Determine orientation of T1 relative to vertical.
                                        
                                        [angle , signv] = orientation(time,m,junctions_array);
                                        angle_before_T1 = [angle_before_T1 angle];
                                        sign_before_T1 = [sign_before_T1 signv];
                                        
                                        
                                        %After T1
                                        %ma = n_of_id(time_after,mid);
                                        %vector_m_after = [junctions_array(time_after,ma).vertex1(1)-junctions_array(time_after,ma).vertex2(1) junctions_array(time_after,ma).vertex1(2)-junctions_array(time_after,ma).vertex2(2) 0];
                                        %angle_test_a(count_t1) = abs(abs(atan2(norm(cross(vector_m_after,vector_vertical)),dot(vector_m_after,vector_vertical))-pi/2)-pi/2);
                                        
                                        
                                        %%%%%%%%%%%%%%%%%%%%%%
                                        
                                    end
                                end
                                
                            end
                            
                            %%%%%%%%%
                            
                        end
                    end
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%test angles before T1 transition

angles_before_T1 = zeros(length(T1_ids),5);

for i = 4:length(T1_ids)
    time = T1_time(i);
    m = n_of_id(time,T1_ids(i));
    for j = 1:5
        t = time + 1 - j;
        if ~isempty(junctions_array(t,m).vertex1) && ~isempty(junctions_array(t,m).vertex2)
            if ~(id_of_n(t,m) == 0)
                [angle] = orientation(t,m,junctions_array);
                angles_before_T1(i,j) = angle;
            end
        end
    end
end

%add last time point before T1
angles_before_T1(:,1) = angle_before_T1';

%only include junctions that weren't at a four-way vertex at least five
%time points back in the determination of angle, by finding nonzero rows
all(angles_before_T1,2)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_size = 25;
%chosen_data = length_array(:,chosen_junction_id);


for time = 1:numberImage %time = window_size+20:numberImage-window_size  %238(95)
    image_junctions = zeros(sizeImage(time,1),sizeImage(time,2));
    %The option 'pic' and 'saveas' are if one wants to create a movie out of
    %the images
    filename = ['image', num2str(time)];
    %pic=figure();
    
    
    for m = 1:numberJunctions(time)
        
        %Colour image of each junction
        for i = 1:length(junctions_array(time,m).junctionCoordinates(:,1))
            image_junctions(junctions_array(time,m).junctionCoordinates(i,1),junctions_array(time,m).junctionCoordinates(i,2),:) = 102;
        end
    end
    
    for chosen_junction_id = stable_T1s([1 2 3 6 7 8 10 11])
        if ~(n_of_id(time,chosen_junction_id) == 0)
            
            %Colour chosen junction
            for i = 1:length(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(:,1))
                image_junctions(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,1),junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,2),:) = 90;
            end
            
            
        end
        
    end
    
    %%%%%%%%%%%%%
    
    
    %         for chosen_junction_id = unique(total_junctions)
    %             if ~(n_of_id(time,chosen_junction_id) == 0)
    %
    %                %Colour chosen junction
    %                for i = 1:length(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(:,1))
    %                    image_junctions(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,1),junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,2),:) = 90;
    %                end
    %             end
    %
    %         end
    
    
    %        count_t = 0;
    %        for chosen_junction_id = T1_ids(find(ismember(T1_time,1:numberImage)))%T1_ids(find(ismember(T1_time,time-20:time+20)))  %shrink_ids(find(ismember(shrink_time,[time-20:time+20])))
    %            count_t = count_t + 1;
    %            T1_transition_time = T1_ids(find(ismember(T1_time,1:numberImage))); %T1_time(find(ismember(T1_time,time-20:time+20)));
    %            T1_transition = T1_transition_time(count_t);
    %            if ~(n_of_id(time,chosen_junction_id) == 0)
    %
    %                %Colour chosen junction
    %                for i = 1:length(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(:,1))
    %                    image_junctions(junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,1),junctions_array(time,n_of_id(time,chosen_junction_id)).junctionCoordinates(i,2),:) = 90;
    %                end
    %
    % %        %Calculate correlation for each of the first neighbours of the T1 junctions
    % %        if ~(n_of_id(time,chosen_junction_id) == 0)
    % %        for n = find(junction_neighbours(time,n_of_id(time,chosen_junction_id),:) == 1)'
    % %            if ~(id_of_n(time,n) == 0)
    % %            chosen_data = length_array(:,chosen_junction_id);
    % %            other_data = length_array(:,id_of_n(time,n));
    % %
    % %            %window the data (colour a junction by crosscorrelation before and
    % %            %after transition)
    % %
    % %            if T1_transition < time
    % %                chosen_data_w = chosen_data(T1_transition-window_size:T1_transition);
    % %                other_data_w = other_data(T1_transition-window_size:T1_transition);
    % %            else
    % %                chosen_data_w = chosen_data(T1_transition:T1_transition+window_size);
    % %                other_data_w = other_data(T1_transition:T1_transition+window_size);
    % %            end
    % %
    % %            [C,lags] = xcorr(chosen_data_w-mean(chosen_data_w),other_data_w-mean(other_data_w),0,'coeff');
    % %            C = 100*(C+1)/2+2;
    % %
    % %
    % %            for i = 1:length(junctions_array(time,n).junctionCoordinates(:,1))
    % %                 image_junctions(junctions_array(time,n).junctionCoordinates(i,1),junctions_array(time,n).junctionCoordinates(i,2),:) = C;
    % %            end
    % %
    % %            end
    % %        end
    % %        end
    %
    %            end
    %
    %        end
    
    
    
    
    %caxis([-101,101])
    cmp = colormap(jet(103)); %hsv
    cmp(1,:) = [0 0 0];
    cmp(103,:) = [1 1 1];
    cmp(102,:) = [0.5 0.5 0.5];
    %image_junctions_sub = image_junctions(170:250,240:320,:);
    %pic = imshow(image_junctions_sub,cmp,'InitialMagnification',200);
    pic = imshow(image_junctions,cmp); %colorbar;
    saveas(pic,filename,'png');
    
    
    
    %     %Plot angle
    %     T1angle = angle_before_T1(4);
    %     angle_sign = sign_before_T1(4);
    %     %The angle_sign is such that for negative sign, the junction is rotated
    %     %in the negative (clockwise) direction relative to the horisontal.
    %     p1 = [junctions_array(time,n_of_id(time,chosen_junction_id)).vertex1(2), junctions_array(time,n_of_id(time,chosen_junction_id)).vertex1(1)];
    %     p2 = [p1(1)-10*cos(T1angle),p1(2)+angle_sign*10*sin(T1angle)]; %[junctions_array(time,n_of_id(time,chosen_junction_id)).vertex2(2), junctions_array(time,n_of_id(time,chosen_junction_id)).vertex2(1)];
    %     line([p1(1) p2(1)],[p1(2) p2(2)],'Color','b','LineWidth',2)
    %     line([500 500],[260 290],'Color','r','LineWidth',2)
    
end