%MSD of vertices
close all
clear all
tic

load('junctions_WT181012.mat')

%%%%%%%%%%%%%% Track vertices %%%%%%

%vertex_array{t} gives the coordinates of all the vertices at time = t
%track_vertices puts the coordinates in a format that the 'track' function
%needs as input
track_vertices = vertex_array{1};
track_vertices(:,3) = 1;

for t = 2:numberImage
    track_vertices = [track_vertices ; vertex_array{t}];
    track_vertices(end-length(vertex_array{t})+1:end,3) = t;
end


dd = length(track_vertices(1,:));
%for mem, generally use 24, but less for Natalia's files
param = struct('mem',2,'dim',dd-1,'good',0,'quiet',0);

vertex_positions = track(track_vertices,12,param);

%Find vertices that persist for the whole movie
vertex_ids = intersect(vertex_positions(find(vertex_positions(:,3) == 1),4),vertex_positions(find(vertex_positions(:,3) == numberImage),4));

%Remove any vertices at the boundary
border1 = 2;
border2 = 2;
border3 = sizeImage(time,1)-1;
border4 = sizeImage(time,2)-1;

border_id = 0;
remove_vertex = 0;
for border_id = 1:length(vertex_ids)
    
    bordervertex = vertex_positions(find(vertex_positions(:,4) == vertex_ids(border_id)),:);
    
    if ~all(bordervertex(:,1) > border1) || ~all(bordervertex(:,1) < border3) || ~all(bordervertex(:,2) > border2) || ~all(bordervertex(:,2) < border4)
        remove_vertex = [remove_vertex; border_id]
    end
end

remove_vertex(1) = [];
vertex_ids(remove_vertex) = [];


%%%%%%%%%%%%Generate vertex coordinates%%%%%%%%%%%%%%%%%%%%%%%
counter = 1;
maxIm = numberImage/2;
vertex_tracks = {};
for v_id = 1:length(vertex_ids)
    walkplot = vertex_positions(find(vertex_positions(:,4) == vertex_ids(v_id)),:);
    
    %only include vertices that move
    if length(unique(walkplot(:,1))) > 1 && length(unique(walkplot(:,2))) > 1
        
        %if length(walkplot(:,2)) > maxIm
        %plot examples of diffusion paths
        
        if counter < 10
            colours = [1:length(walkplot(:,1))]';
            z = zeros(size(walkplot(:,1)));
            subplot(3,3,counter), surface([walkplot(:,1),walkplot(:,1)],[walkplot(:,2),walkplot(:,2)],[z,z],[colours,colours],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
            pause(1.0)
        end
        
        vx = walkplot(:,1);
        vy = walkplot(:,2);
        
        vertex_tracks{counter} = [[1:length(vx)]'*30  vx*0.08961  vy*0.08961];
        
        counter = counter + 1;
        %end
    end
    
end


%%%%%%%%%%Calculate MSD%%%%%%%%%%%%%%%

analyze_MSD(vertex_tracks(19))
%analyze_msd_directed_motion(vertex_tracks)

% dt = 30;
% figure()
% nval= [1:1:time2-1];
% plot(dt*nval,(0.08961^2)*msd_plot)
% hold on
% a = polyfit(dt*nval(1:end-20),(0.08961^2)*msd_plot(1:end-20),1)
% plot(dt*nval,a(1)*dt*nval+a(2),'g')
%
% figure()
% loglog(dt*nval,(0.08961^2)*msd_plot)
% hold on
% loglog(dt*nval,dt*nval,'g')
toc

