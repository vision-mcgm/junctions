close all
clear all

load('junction_data.mat')


range = 3;
   
%%%%%%%%%%%%%%%%Before T1 transition
count_t = 1;
count_t1 = 1;
for chosen_junction_id = T1_ids 
    T1_transition = T1_time(count_t);
    count_t = count_t + 1;
    
    %Before transition    
    time = T1_transition - range;
    if time > 0
    if ~(n_of_id(time,chosen_junction_id) == 0)
    
                m = n_of_id(time,chosen_junction_id);
                vector_m = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
                vector_vertical = [0 1 0];
        
                %choose vector m to point upwards (take vertex with largest y-value
                %and substract vertex with lowest y-value from it to form the
                %vector)
                max_vertex = find(max([junctions_array(time,m).vertex1(2) junctions_array(time,m).vertex2(2)]));
                vectors(1,:) = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
                vectors(2,:) = [junctions_array(time,m).vertex2(1)-junctions_array(time,m).vertex1(1) junctions_array(time,m).vertex2(2)-junctions_array(time,m).vertex1(2) 0];
                vector_m = vectors(max_vertex,:);

                %find angle vector makes with vertical (between 0 an pi/2)
                angle_vertical = abs(abs(atan2(norm(cross(vector_m,vector_vertical)),dot(vector_m,vector_vertical))-pi/2)-pi/2); 
        
                %get sign of angle by calculating the crossproduct between
                %vector_m and vector_vertical - this gives a normal to the plane.
                %If if points upwards (positive y), vector m is a clockwise
                %rotation relative to vector_vertical and the angle theta by which
                %to rotate the plane is positive (counterclockwise)
                v_normal = cross(vector_m,vector_vertical);
                angle_vertical = sign(v_normal(3))*angle_vertical;
                
                angle_before_T1(count_t1) = angle_vertical;
                
                count_t1 = count_t1 + 1;
    end
    end
end
     
range = 3;
   
%%%%%%%%%%%%%%%%After T1 transition
count_t = 1;
count_t1 = 1;
for chosen_junction_id = T1_ids 
    T1_transition = T1_time(count_t);
    count_t = count_t + 1;
    
    %Before transition    
    time = T1_transition + range;
    if time < numberImage
    if ~(n_of_id(time,chosen_junction_id) == 0)
    
                m = n_of_id(time,chosen_junction_id);
                vector_m = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
                vector_vertical = [0 1 0];
        
                %choose vector m to point upwards (take vertex with largest y-value
                %and substract vertex with lowest y-value from it to form the
                %vector)
                max_vertex = find(max([junctions_array(time,m).vertex1(2) junctions_array(time,m).vertex2(2)]));
                vectors(1,:) = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
                vectors(2,:) = [junctions_array(time,m).vertex2(1)-junctions_array(time,m).vertex1(1) junctions_array(time,m).vertex2(2)-junctions_array(time,m).vertex1(2) 0];
                vector_m = vectors(max_vertex,:);

                %find angle vector makes with vertical (between 0 an pi/2)
                angle_vertical = abs(abs(atan2(norm(cross(vector_m,vector_vertical)),dot(vector_m,vector_vertical))-pi/2)-pi/2); 
        
                %get sign of angle by calculating the crossproduct between
                %vector_m and vector_vertical - this gives a normal to the plane.
                %If if points upwards (positive y), vector m is a clockwise
                %rotation relative to vector_vertical and the angle theta by which
                %to rotate the plane is positive (counterclockwise)
                v_normal = cross(vector_m,vector_vertical);
                angle_vertical = sign(v_normal(3))*angle_vertical;
                
                angle_after_T1(count_t1) = angle_vertical;
                
                count_t1 = count_t1 + 1;
    end
    end
end

