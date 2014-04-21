%Get data on position and angle for each junction, along with the
%correlations and neighbour information
%saved as 'angle_correlation.mat'

clear all
close all

load('junction_data.mat')

%Parameter. Width of window over which crosscorrelation is calculated
window_size = 25;

%Preallocate variables
crosscorrelation_data = NaN(numberImage,max(numberJunctions),max(numberJunctions));
coordinate_data_x = NaN(numberImage,max(numberJunctions),max(numberJunctions));
coordinate_data_y = NaN(numberImage,max(numberJunctions),max(numberJunctions));
angle_data = NaN(numberImage,max(numberJunctions),max(numberJunctions));
distance_data = NaN(numberImage,max(numberJunctions),max(numberJunctions));
angle_v_vector = NaN(1,max(numberJunctions));


time_counter = 1;
%Use non-overlapping windows. Only sample time points separated by the window size.
for time = 1:window_size:numberImage-window_size
    
    display(time)
    
    %For each junction, get the crosscorrelation values for all other
    %junctions along with their angle and distance relative to that junction
    
    for m = 1:numberJunctions(time)
    
    %only do this for junctions with two vertices at both start and end of time window to exclude border
    %junctions
    if ~(id_of_n(time,m) == 0) && ~(n_of_id(time + window_size,id_of_n(time,m)) == 0)
        mid = id_of_n(time,m);
    if ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2) && ~isempty(junctions_array(time+window_size,n_of_id(time+window_size,mid)).vertex1) && ~isempty(junctions_array(time+window_size,n_of_id(time+window_size,mid)).vertex2)
    
        %get crosscorrelation for junction m relative to all other
        %junctions in image
        for n = 1:numberJunctions(time)
        
        %exclude calculation of junction with itself
        if ~(m == n)
         
        %only get crosscorrelation for junctions that also have two
        %vertices      
        if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
            nid = id_of_n(time,n);
        if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) && ~isempty(junctions_array(time+window_size,n_of_id(time+window_size,nid)).vertex1) && ~isempty(junctions_array(time+window_size,n_of_id(time+window_size,nid)).vertex2)  

        data1 = length_array(:,mid);
        data2 = length_array(:,nid);
        
        
        %window the data (should the time point be in the middle of the
        %window instead?)
        data1_w = data1(time:time+window_size);
        data2_w = data2(time:time+window_size);
        
        [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
        crosscorrelation_data(time_counter,m,n) = C;
    
        %To calculate angle, create 3D vectors.
        %vector_m = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
        %Calculate angle vector_m makes with vertical
        vector_vertical = [0 1 0];
        
        %Choose vector_m to point upwards (take vertex with largest y-value
        %and substract vertex with lowest y-value from it to form the
        %vector)
        max_vertex = find(max([junctions_array(time,m).vertex1(2) junctions_array(time,m).vertex2(2)]));
        vectors(1,:) = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
        vectors(2,:) = [junctions_array(time,m).vertex2(1)-junctions_array(time,m).vertex1(1) junctions_array(time,m).vertex2(2)-junctions_array(time,m).vertex1(2) 0];
        vector_m = vectors(max_vertex,:);

        %Find smallest angle vector makes with vertical (between 0 an pi/2)
        angle_vertical = abs(abs(atan2(norm(cross(vector_m,vector_vertical)),dot(vector_m,vector_vertical))-pi/2)-pi/2); 
        
        %get sign of angle by calculating the crossproduct between
        %vector_m and vector_vertical - this gives a normal to the plane.
        %If the normal points upwards (positive z), vector m is a clockwise
        %rotation relative to vector_vertical and the angle theta by which
        %to rotate the plane is positive (counterclockwise)
        v_normal = cross(vector_m,vector_vertical);
        angle_vertical = sign(v_normal(3))*angle_vertical;
        
        %Define vector for other junction
        vector_n = [junctions_array(time,n).vertex1(1)-junctions_array(time,n).vertex2(1) junctions_array(time,n).vertex1(2)-junctions_array(time,n).vertex2(2) 0];
        %Calculate smallest angle between vector_m and vector_n (between 0
        %and pi/2)
        angle_mn = abs(abs(atan2(norm(cross(vector_m,vector_n)),dot(vector_m,vector_n))-pi/2)-pi/2); 
        %Get euclidian distance between midpoints of junctions
        distance = sqrt((junctions_array(time,m).midpoint(1)-junctions_array(time,n).midpoint(1))^2+(junctions_array(time,m).midpoint(2)-junctions_array(time,n).midpoint(2))^2);
        
        %Rotate whole coordinate system such that junction m is aligned
        %vertically and has midpoint at (0,0). Rotate the xy-plane
        %by the angle theta about the z-axis. The coordinates of a point are
        %transformed as x' = cos(theta)x+sin(theta)y, y' =
        %-sin(theta)x+cos(theta)y. The angle theta is the angle the
        %junction m makes with the vertical, as this will shift the entire
        %coordinate system such that the junction m is aligned vertically.
        %The transformation is applied to the midpoint of vertex n

        %Define point to rotate
        v_before = [junctions_array(time,m).midpoint(1) junctions_array(time,n).midpoint(1); junctions_array(time,m).midpoint(2) junctions_array(time,n).midpoint(2)];
        
        center = repmat([junctions_array(time,m).midpoint(1); junctions_array(time,m).midpoint(2)], 1, 2);
        %Define rotation matrix
        theta = angle_vertical;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        %Do the rotation - first shift points in the plane so midpoint
        %(center of rotation) is at the origin
        v_compare = v_before - center;
        v_trans = R*(v_before - center);
        %Plot midpoints before and after rotation to check rotation is
        %correct
%         figure()
%         display(angle_vertical*180/pi)
%         plot([0 v_compare(1,2)],[0 v_compare(2,2)],'g')
%         hold on
%         plot([0 v_trans(1,2)],[0 v_trans(2,2)],'b')
%         pause(5)
%         close all
        
        %Store transformed midpoints of n
        coordinate_data_x(time_counter,m,n) = [v_trans(1,2)];
        coordinate_data_y(time_counter,m,n) = [v_trans(2,2)];

        angle_data(time_counter,m,n) = angle_mn;
        distance_data(time_counter,m,n) = distance;
        
        end
        end
        end
        end
    end
    end
    end
    time_counter = time_counter + 1;
end

%Average and plot the results
load('angle_distance_data.mat')

%Create 2D array of bins to put crosscorrelation values into

n_bins = 100;

%get limits of plot
min_x = min(coordinate_data_x(:));
min_y = min(coordinate_data_y(:));

max_x = max(coordinate_data_x(:));
max_y = max(coordinate_data_y(:));

length_x = abs(max_x)+abs(min_x);
length_y = abs(max_y)+abs(min_y);

%Define bins such that the point (0,0) is at a bin edge (this solution
%means that the last (largest value) bin in the x and y directions can vary
%in size. Since we do not take the outermost bins into account, this is not
%a problem
bin_width_x = abs(min_x)/(n_bins/2);
bin_width_y = abs(min_y/(n_bins/2));
bin_midpoints_x = [min_x+bin_width_x/2:bin_width_x:max_x-bin_width_x/2];
bin_midpoints_y = [min_y+bin_width_y/2:bin_width_y:max_y-bin_width_y/2];
bin_edges_x = [min_x:bin_width_x:max_x];
bin_edges_y = [min_y:bin_width_y:max_y];

%Preallocate variables
bins = NaN(n_bins,n_bins,max(numberJunctions));
average_bins= zeros(n_bins,n_bins);
sum_bins = zeros(n_bins,n_bins);
number_in_bins = zeros(n_bins,n_bins);

%Sort correlation_values into bins
time_counter = 1;

for time = 1:window_size:numberImage-window_size;
    
    display(time)
    
for m = 1:numberJunctions(time)
      
    counter = 1;
    bins = NaN(n_bins,n_bins,max(numberJunctions));
  
for n = 1:numberJunctions(time)
    
    %exclude cases of m = n (otherwise you get an strong correlation in the
    %middle
    
    if ~(m == n)
    
    bin_data_x = zeros(length(bin_edges_x)+1);
    bin_data_y = zeros(length(bin_edges_y)+1);
    index1 = 0;
    index2 = 0;

    if ~(isnan(coordinate_data_x(time_counter,m,n)))
        
        %Sort the data into bins. Add the data coordinate to the vector
        %bin_edges, then sort the vector and find the where the data
        %coordinate is ranked - this will be the bin the data point should
        %be in
        bin_data_x = [bin_edges_x coordinate_data_x(time_counter,m,n)];
 
        index1 = max(find(sort(bin_data_x) == coordinate_data_x(time_counter,m,n))-1);
        
        bin_data_y = [bin_edges_y coordinate_data_y(time_counter,m,n)];
 
        index2 = max(find(sort(bin_data_y) == coordinate_data_y(time_counter,m,n))-1);      
        
        bins(index1,index2,counter) = crosscorrelation_data(time_counter,m,n);
        counter = counter + 1;
    end
    end
    
end

for i1 = 1:n_bins
    for i2 = 1:n_bins
        sum_bins(i1,i2) = sum_bins(i1,i2) + nansum(bins(i1,i2,:));
        number_in_bins(i1,i2) = number_in_bins(i1,i2) + length(find(~(isnan(bins(i1,i2,:)))));
    end
end

end

time_counter = time_counter+1;

end

for i1 = 1:n_bins
    for i2 = 1:n_bins
        average_bins(i1,i2) = sum_bins(i1,i2)/number_in_bins(i1,i2);
    end
end

%load('heatmap.mat')

figure()
surf(bin_midpoints_y(1:99),bin_midpoints_x(1:99),average_bins(1:99,1:99),'EdgeColor','none')
caxis([-0.38 0.38])
shading flat %interp
hold on
plot([0 0],[52/2 -52/2],'k','LineWidth',2)
figure()
surf(bin_midpoints_y(41:61),bin_midpoints_x(41:61),average_bins(41:61,41:61),'EdgeColor','none')
caxis([-0.3 0.3])
shading flat
hold on
plot([0 0],[52/2 -52/2],'k')


%Problemer: Hvorfor er plottet ikke symmetrisk (selvom hvis man bruger interpolate option for shading bliver det) omkring y-aksen? Hvorfor er
%det orienteret modsat af hvad det burde vaere. Hvorfor er effekten saa
%lille at hvis det bliver plotter kun for et tidskridt, er der naesten
%ingen effekt?
