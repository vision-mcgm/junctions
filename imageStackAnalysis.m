clear all
close all
clc

tic

numberJunctionsMany = 1000;
numberImage = 10; %299; %82; %(317)change this to number of images in folder
numberJunctions = zeros(numberImage,1);
sizeImage = zeros(numberImage,3);

%Initialize array of objects - For objects the pre-allocation works by assigning one of the objects to the very last field in the array. Matlab then fills the other fields before that with objects (handles) that it creates by calling the constructor of that object with no arguments.
junctions_array(numberImage,numberJunctionsMany) = junction();
for t = 1:numberImage
    for n = 1:numberJunctionsMany
        junctions_array(t,n) = junction();
    end
end

%Initialize array of vertex_arrays
vertex_array = cell(numberImage,1);

ffolder = '181012 SPE2 bd_fate files'; %'GBE_OFA_stable';%'bdfate GBE 240713 z4';%'181012 SPE2 bd_fate files';%'bd_fate_high_res'; %'DECADGFP ROKiR 240113 SPE2 bd_fate'; %'bd_fate_T1' to use other folder
fname = 'bd_fate%.3d.png';

for time = 1:numberImage %numberImage
    
    display(time)
    
    %clear workspace for the temporary variables and the global variables
    %used in the function for the previous image
    clear junctions_holder vertex_holder numJunctions_holder;
    clearvars -global;
    
    %time_t = 1000 + time; %since the images are named from 1000
    
    filename = fullfile(ffolder,fname);
    %im = imread(sprintf(filename,time));
    im = imread([filename(1:end-8) num2str(time,'%03d') '.png']); %NB The other line does not work on Jeppes computer. This does.
    
    sizeImage(time,:) = size(im);
    
    [junctions_holder, vertex_holder, numJunctions_holder] = findEdges(im);
    vertex_array{time} = vertex_holder;
    numberJunctions(time) = numJunctions_holder;
    %remove empty entries
    for n = 1:numberJunctions(time)
        dimension = size(junctions_holder(n).junctionCoordinates);
        if dimension(1) < 1
            junctions_holder(n) = [];
        end
    end
    %Assign junctions to array
    junctions_array(time,1:length(junctions_holder)) = junctions_holder;
    
end
toc

%Truncate junctions_array
junctions_array = junctions_array(:,1:max(numberJunctions));

%Only save junctions_array, vertex_array, numberJunctions, numberImage
savename = [ffolder '_junctions' '.mat'];
save(savename,'junctions_array','vertex_array','numberJunctions','numberImage','sizeImage')

