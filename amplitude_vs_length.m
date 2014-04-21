close all 
clear all

load('junctions_rok.mat')

%Exclude border junctions, by requiring junctions to have at least two
%vertices


%Take variance of time series (take different lengths of same time series
%to demonstrate that variance depends on the length of the time series,
%then do variance for the same length of different time series)

time1 = 2;
time2 = 420;
count = 1;

for n = 1:numberJunctions(time1)
if ~isempty(junctions_array(time1,n).vertex1) && ~isempty(junctions_array(time1,n).vertex2) && ~isempty(junctions_array(time2,n).vertex1) && ~isempty(junctions_array(time2,n).vertex2)
if ~(id_of_n(time1,n) == 0) 
    if ~(n_of_id(time2,id_of_n(time1,n)) == 0)
    
    display (id_of_n(time1,n))
    
    idj = id_of_n(time1,n);
    
    ts = length_array(time1:time2,idj);
    ts_m = ts - mean(ts);
    
    %Pick length of time window such that data is stationary (0 indicates that
    %you cannot reject null hypothesis that the time series is unit root).
    stat_test(count) = adftest(ts_m);
    
    ts_mean(count) = mean(ts);
    ts_var(count) = var(ts_m);
    id_of_count(count) = idj;
    count = count+1;
   
    end

end
end
end
    
figure()
scatter(ts_mean,ts_var)
%ylim([0 200])
xlabel('Mean junction length')
ylabel('Variance')


test0 = find(stat_test == 0)
figure()
scatter(ts_mean(test0),ts_var(test0))
%ylim([0 200])

test1 = find(stat_test == 1)
figure()
scatter(ts_mean(test1),ts_var(test1))
%ylim([0 200])


figure()
hist(stat_test(:))


    
    



    if ~(n_of_id(time,chosen_junction_id) == 0)
        for n = find(junction_neighbours(time,n_of_id(time,chosen_junction_id),:) == 1)'
            
                chosen_data = length_array(:,chosen_junction_id);
                other_data = length_array(:,id_of_n(time,n)); 
                m = n_of_id(time,chosen_junction_id);
                



for range = 1:25

count_before = 1;    
%%%%%%%%%%%%%%%%Before T1 transition
count_t = 1;
for chosen_junction_id = T1_ids 
    T1_transition = T1_time(count_t);
    count_t = count_t + 1;
    
    if T1_transition > window_size
    %Before transition    
    time = T1_transition - range;
    %Calculate correlation for each of the first neighbours of the T1
    if ~(n_of_id(time,chosen_junction_id) == 0)
        for n = find(junction_neighbours(time,n_of_id(time,chosen_junction_id),:) == 1)'
            if ~(id_of_n(time,n) == 0)
                chosen_data = length_array(:,chosen_junction_id);
                other_data = length_array(:,id_of_n(time,n)); 
                m = n_of_id(time,chosen_junction_id);
                
                if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                
                chosen_data_w = chosen_data(T1_transition-window_size:T1_transition);      
                other_data_w = other_data(T1_transition-window_size:T1_transition); 
                [C,lags] = xcorr(chosen_data_w-mean(chosen_data_w),other_data_w-mean(other_data_w),0,'coeff');        
                crosscorrelation_before(range,count_before) = C;
    
                %Calculate angle. Approximate each junction to be a
                %straight line between the two vertices. For neighbouring
                %junctions use cosine rule. For all other junctions, find
                %the smallest angle between them (will be between 0 and
                %pi/2)
      
                
                vector_m = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
                vector_n = [junctions_array(time,n).vertex1(1)-junctions_array(time,n).vertex2(1) junctions_array(time,n).vertex1(2)-junctions_array(time,n).vertex2(2) 0];            
                                              
                %Using the cosine rule, calculate the angle corresponding
                %to the vertex shared by the two junctions
                vertices = [junctions_array(time,m).vertex1 ; junctions_array(time,m).vertex2 ; junctions_array(time,n).vertex1 ; junctions_array(time,n).vertex2];
                %find shared vertex
                
                    [sorted_vertices,~,uniqueID] = unique(vertices,'rows');
                    p_shared = sorted_vertices(mode(uniqueID),:);
                    %remove the index of the shared vertex from the uniqueID
                    %list
                    uniqueID(uniqueID == mode(uniqueID)) = [];
                    %define the other two vertices
                    p_1 = sorted_vertices(uniqueID(1),:);
                    p_2 = sorted_vertices(uniqueID(2),:);
                    %define line segments from the shared vertex to each of the
                    %other vertices
                    line_seg1 = [p_1(1)-p_shared(1)  p_1(2)-p_shared(2) 0];
                    line_seg2 = [p_2(1)-p_shared(1) p_2(2)-p_shared(2) 0];
                    %calculate angle
                    angle_before(range,count_before) = acos(dot(line_seg1/norm(line_seg1),line_seg2/norm(line_seg2)));                   
                    count_before = count_before + 1;
                end
                end
 
            end
        end
    end 

end
end










%Plot variance as a function of mean length of junction