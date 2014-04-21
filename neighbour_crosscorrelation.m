close all
clear all

load('junction_data.mat')

window_size = 25;

crosscorrelation_neighbours = NaN(max(numberJunctions),1); 

first_neighbours_c = [];
second_neighbours_c = [];
third_neighbours_c = [];
fourth_neighbours_c = [];
fifth_neighbours_c = [];
            

for time = 1:window_size:numberImage-window_size
    display(time)
    for m = 1:numberJunctions(time)

        first_neighbours = find(junction_neighbours(time,m,:) == 1);
        
        second_neighbours = [];
        for k = 1:length(first_neighbours)
        s_neighbours = find(junction_neighbours(time,first_neighbours(k),:) == 1);
        second_neighbours = [second_neighbours; s_neighbours];
        end

        %remove both the junction m and any junctions that are also first
        %neighbours
        second_neighbours(second_neighbours == m) = [];
        second_neighbours(find(ismember(second_neighbours,first_neighbours) == 1)) = [];
       
        
        third_neighbours = [];
        for k = 1:length(second_neighbours)
        t_neighbours = find(junction_neighbours(time,second_neighbours(k),:) == 1);
        third_neighbours = [third_neighbours; t_neighbours];
        end

        third_neighbours(third_neighbours == m) = [];
        third_neighbours(find(ismember(third_neighbours,first_neighbours) == 1)) = [];
        third_neighbours(find(ismember(third_neighbours,second_neighbours) == 1)) = [];

        
        fourth_neighbours = [];
        for k = 1:length(third_neighbours)
        f_neighbours = find(junction_neighbours(time,third_neighbours(k),:) == 1);
        fourth_neighbours = [fourth_neighbours; f_neighbours];
        end

        fourth_neighbours(fourth_neighbours == m) = [];
        fourth_neighbours(find(ismember(fourth_neighbours,first_neighbours) == 1)) = [];
        fourth_neighbours(find(ismember(fourth_neighbours,second_neighbours) == 1)) = [];
        fourth_neighbours(find(ismember(fourth_neighbours,third_neighbours) == 1)) = [];
            
        
        
        fifth_neighbours = [];
        for k = 1:length(fourth_neighbours)
        fi_neighbours = find(junction_neighbours(time,fourth_neighbours(k),:) == 1);
        fifth_neighbours = [fifth_neighbours; fi_neighbours];
        end

        fifth_neighbours(fifth_neighbours == m) = [];
        fifth_neighbours(find(ismember(fifth_neighbours,first_neighbours) == 1)) = [];
        fifth_neighbours(find(ismember(fifth_neighbours,second_neighbours) == 1)) = [];
        fifth_neighbours(find(ismember(fifth_neighbours,third_neighbours) == 1)) = [];
        fifth_neighbours(find(ismember(fifth_neighbours,fourth_neighbours) == 1)) = [];
            
        
        sixth_neighbours = [];
        for k = 1:length(fifth_neighbours)
        si_neighbours = find(junction_neighbours(time,fifth_neighbours(k),:) == 1);
        sixth_neighbours = [sixth_neighbours; si_neighbours];
        end

        sixth_neighbours(sixth_neighbours == m) = [];
        sixth_neighbours(find(ismember(sixth_neighbours,first_neighbours) == 1)) = [];
        sixth_neighbours(find(ismember(sixth_neighbours,second_neighbours) == 1)) = [];
        sixth_neighbours(find(ismember(sixth_neighbours,third_neighbours) == 1)) = [];
        sixth_neighbours(find(ismember(sixth_neighbours,fourth_neighbours) == 1)) = [];
        sixth_neighbours(find(ismember(sixth_neighbours,fifth_neighbours) == 1)) = [];
        
        seventh_neighbours = [];
        for k = 1:length(sixth_neighbours)
        se_neighbours = find(junction_neighbours(time,sixth_neighbours(k),:) == 1);
        seventh_neighbours = [seventh_neighbours; se_neighbours];
        end

        seventh_neighbours(seventh_neighbours == m) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,first_neighbours) == 1)) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,second_neighbours) == 1)) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,third_neighbours) == 1)) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,fourth_neighbours) == 1)) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,fifth_neighbours) == 1)) = [];
        seventh_neighbours(find(ismember(seventh_neighbours,sixth_neighbours) == 1)) = [];
            
        
        eigth_neighbours = [];
        for k = 1:length(seventh_neighbours)
        ei_neighbours = find(junction_neighbours(time,seventh_neighbours(k),:) == 1);
        eigth_neighbours = [eigth_neighbours; ei_neighbours];
        end

        eigth_neighbours(eigth_neighbours == m) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,first_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,second_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,third_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,fourth_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,fifth_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,sixth_neighbours) == 1)) = [];
        eigth_neighbours(find(ismember(eigth_neighbours,seventh_neighbours) == 1)) = [];
        
        
    %calculate crosscorrelation           
    if ~isempty(first_neighbours)
    crosscorrelation_neighbours = NaN(max(numberJunctions),1);        
    count = 1;
    for n = first_neighbours'
        
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                               

            data1 = length_array(:,id_of_n(time,n));
            data2 = length_array(:,id_of_n(time,m));
        
        
            %window the data
            data1_w = data1(time:time+window_size);
            data2_w = data2(time:time+window_size);
        
            [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
            crosscorrelation_neighbours(count) = C;            
            
            count = count + 1;
            
            end
            end
                  
    end
            if ~(length(first_neighbours) == 0)
            first_neighbours_c = [first_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            end 
    
    end
    
    if ~isempty(second_neighbours)
    crosscorrelation_neighbours = NaN(max(numberJunctions),1); 
    count = 1;
    for n = second_neighbours'
        
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                               

            data1 = length_array(:,id_of_n(time,n));
            data2 = length_array(:,id_of_n(time,m));
        
        
            %window the data
            data1_w = data1(time:time+window_size);
            data2_w = data2(time:time+window_size);
        
            [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
            crosscorrelation_neighbours(count) = C;            
            
            count = count + 1;
            
            end
            end         
    end
            if ~(length(second_neighbours) == 0)
            second_neighbours_c = [second_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            end 
    end
    
    if ~isempty(third_neighbours)
    crosscorrelation_neighbours = NaN(max(numberJunctions),1); 
    count = 1;
    for n = third_neighbours'
        
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                               

            data1 = length_array(:,id_of_n(time,n));
            data2 = length_array(:,id_of_n(time,m));
        
        
            %window the data
            data1_w = data1(time:time+window_size);
            data2_w = data2(time:time+window_size);
        
            [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
            crosscorrelation_neighbours(count) = C;            
            
            count = count + 1;
            
            end
            end
    end
            if ~(length(third_neighbours) == 0)
            third_neighbours_c = [third_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            end  
    end
    
    if ~isempty(fourth_neighbours)
    crosscorrelation_neighbours = NaN(max(numberJunctions),1); 
    count = 1;
    for n = fourth_neighbours'
        
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                               

            data1 = length_array(:,id_of_n(time,n));
            data2 = length_array(:,id_of_n(time,m));
        
        
            %window the data
            data1_w = data1(time:time+window_size);
            data2_w = data2(time:time+window_size);
        
            [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
            crosscorrelation_neighbours(count) = C;            
            
            count = count + 1;
            
            end
            end
                   
    end
    
            if ~(length(fourth_neighbours) == 0)
            fourth_neighbours_c = [fourth_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            end  
    end
    
    if ~isempty(fifth_neighbours)
    crosscorrelation_neighbours = NaN(max(numberJunctions),1); 
    count = 1;
    for n = fifth_neighbours'
        
            %only carry out crosscorrelation and angle calculation is both junctions n and m have two vertices    
            if ~isempty(junctions_array(time,n).vertex1) && ~isempty(junctions_array(time,n).vertex2) &&  ~isempty(junctions_array(time,m).vertex1) && ~isempty(junctions_array(time,m).vertex2)
                %only carry out crosscorrelation and angle calculation if
                %both junctions still exist at the end of the time window
                %for which the crosscorrelation is carried out
                if ~(id_of_n(time,n) == 0) && ~(id_of_n(time,m) == 0) && ~(n_of_id(time+window_size,id_of_n(time,n)) == 0) && ~(n_of_id(time+window_size,id_of_n(time,m)) == 0)
                               

            data1 = length_array(:,id_of_n(time,n));
            data2 = length_array(:,id_of_n(time,m));
        
        
            %window the data
            data1_w = data1(time:time+window_size);
            data2_w = data2(time:time+window_size);
        
            [C,lags] = xcorr(data1_w-mean(data1_w),data2_w-mean(data2_w),0,'coeff');
        
            crosscorrelation_neighbours(count) = C;            
            
            count = count + 1;
            
            end
            end
                     
    end 
            if ~(length(fifth_neighbours) == 0)
            fifth_neighbours_c = [fifth_neighbours_c ; crosscorrelation_neighbours(1:count-1)];
            end 
    
    end
    
    end
    
end


%%%%%%%%%%%%plot histogram

subplot(2,2,1), hist(first_neighbours_c,100)
title('Crosscorrelation first neighbours')
subplot(2,2,2), hist(second_neighbours_c,100)
title('Crosscorrelation second neighbours')
subplot(2,2,3), hist(third_neighbours_c,100)
title('Crosscorrelation third neighbours')
subplot(2,2,4), hist(fourth_neighbours_c,100)
title('Crosscorrelation fourth neighbours')



%%%%%%%%%%%%make boxplot
n1 = first_neighbours_c;
n2 = second_neighbours_c;
n3 = third_neighbours_c;
n4 = fourth_neighbours_c;
n5 = fifth_neighbours_c;
%n6 = sixth_neighbours_c';
%n7 = seventh_neighbours_c';
%n8 = eigth_neighbours_c';


group = [repmat({'1. neighbours'}, length(n1), 1); repmat({'2. neighbours'}, length(n2), 1); repmat({'3. neighbours'}, length(n3), 1); repmat({'4. neighbours'}, length(n4), 1); repmat({'5. neighbours'}, length(n5), 1) ];
boxplot([n1;n2;n3;n4;n5], group)



