close all
clear all

load('junctions_highres.mat')

%Find the crosscorrelation and angles for a junction and it's first,
%second, and third neighbours

%The length of the time window over which the correlation is calculated
window_size = 25;

crosscorrelation_neighbours = NaN(max(numberJunctions),1); 

first_neighbours_c = [];
angle_first_neighbours = [];


for time = 1:window_size:numberImage-window_size
    display(time)
    
    for m = 1:numberJunctions(time)
                
        first_neighbours = find(junction_neighbours(time,m,:) == 1);
             
            %calculate crosscorrelation
            if ~isempty(first_neighbours)
            count = 1;
            for n = first_neighbours' %; second_neighbours; third_neighbours]'
            
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                 

                data1 = length_array(:,id_of_n(time,n));
                data2 = length_array(:,id_of_n(time,m));
        
        
                %window the data - current time is start point for window
                %(this could be changed to be midpoint)
                data1_w = data1(time:time+window_size);
                data2_w = data2(time:time+window_size);
                
                %Calculating the crosscorrelation for zero lag
                [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
                crosscorrelation_neighbours(count) = C;            
            
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
                    angle(count) = acos(dot(line_seg1/norm(line_seg1),line_seg2/norm(line_seg2)));                   
                
                count = count + 1;
            
                end
            end
            end
            
            if ~(length(first_neighbours) == 0)
            first_neighbours_c = [first_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            angle_first_neighbours = [angle_first_neighbours angle(1:count-1)];
            end
            
            end
                              
          
    end
    
end


figure()
hist(angle_first_neighbours*180/pi,100)


% for v = 1:length(angle_first_neighbours)
% Lecuit_angle(v) = acos((0.5-cos(angle_first_neighbours(v)))/((5/4-cos(angle_first_neighbours(v)))^(1/2)));
% end
% figure()
% [f,x] = hist(2*Lecuit_angle*180/pi,[30,45,60,75,90,105,120,135,150])
% bar(x,f/sum(f))
% xlim([37.5 127.5])

%Sort correlation values into bins according to their angle

n_bins = 100; 
bins = NaN(n_bins,length(first_neighbours_c));
average_bins= NaN(n_bins,1);

length_a = pi;

bin_width = length_a/n_bins;
bin_edges = [0:bin_width:length_a-bin_width];
bin_midpoints = [0+bin_width/2:bin_width:length_a-bin_width/2];

%sort correlation_values into bins

counter = 1;
counting = ones(n_bins,1);
sum_bins = zeros(n_bins,1);
number_in_bins = zeros(n_bins,1);

for k = 1:length(first_neighbours_c)
    
    bin_data_a = [];
    index_angle = 0;
    
    bin_data_a = [bin_edges angle_first_neighbours(counter)];
 
    index_angle = max(find(sort(bin_data_a) == angle_first_neighbours(counter))-1);
  
    sum_bins(index_angle) = sum_bins(index_angle)+first_neighbours_c(counter); 
    
    for i1 = index_angle
        bins(i1,counting(i1)) = first_neighbours_c(counter);
        number_in_bins(i1) = number_in_bins(i1) + 1;
        counting(i1) = counting(i1) + 1;
    end
    
    counter = counter + 1;    
end


for i1 = 1:n_bins
average_bin(i1) = sum(sum_bins(i1))/sum(number_in_bins(i1));
standard_deviation(i1) = nanstd(bins(i1,:))/sqrt(number_in_bins(i1));
end

%plot only the interval for which there are more than 50 datapoints in the
%bins
figure()
%plot(bin_midpoints(42:88)*180/pi,average_bin(42:88))
errorbar(bin_midpoints(45:87)*180/pi,average_bin(45:87),standard_deviation(45:87),'xb')
%axis([80,160,-0.5,0.2])
xlabel('Angle (degrees)')
ylabel('Crosscorrelation value, lag 0')

figure()
plot(bin_midpoints*180/pi,number_in_bins)
