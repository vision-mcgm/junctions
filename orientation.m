function [ angle , signv] = orientation(time,m,junctions_array)
%Find orientation of junction relative to horisontal (as viewed when
%plotting using imshow, but note that this has the x- and y-coordinates
%swapped, so in cartesians it is relative to vertical).

    vector_m = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
    vector_vertical = [0 1 0];
   
                               
    %choose vector m to point upwards (take vertex with largest y-value
    %and substract vertex with lowest y-value from it to form the vector)
    max_vertex = find(max([junctions_array(time,m).vertex1(2) junctions_array(time,m).vertex2(2)]));
    vectors(1,:) = [junctions_array(time,m).vertex1(1)-junctions_array(time,m).vertex2(1) junctions_array(time,m).vertex1(2)-junctions_array(time,m).vertex2(2) 0];
    vectors(2,:) = [junctions_array(time,m).vertex2(1)-junctions_array(time,m).vertex1(1) junctions_array(time,m).vertex2(2)-junctions_array(time,m).vertex1(2) 0];
    vector_m = vectors(max_vertex,:);

    %find angle vector makes with vertical (between 0 an pi/2)
    angle_vertical = abs(abs(atan2(norm(cross(vector_m,vector_vertical)),dot(vector_m,vector_vertical))-pi/2)-pi/2); 
    angle = angle_vertical;   
    
    %get sign of angle by calculating the crossproduct between
    %vector_m and vector_vertical - this gives a normal to the plane.
    %If it points upwards (positive y), vector m is a clockwise
    %rotation relative to vector_vertical and the angle theta by which
    %to rotate the plane is positive (counterclockwise)
    v_normal = cross(vector_m,vector_vertical);
    signv = sign(v_normal(3)); %*angle_vertical;

end

%could just use tan^-1!