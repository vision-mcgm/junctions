%Create histograms of the frequency of T1 and potential T1 events.
%Create histograms of the time interval between successive events.

%load('')

%%%%Distributions of events per junction and time between events

%Total number of junctions included in analysis is

totalj = length(unique(total_junctions));

%Histogram of the number of four-way vertex contractions per junction
almostT1count = 0;
[almostT1count,un]=hist(almost_T1_ids,unique(almost_T1_ids));

figure()
hist(almostT1count) %include the junctions with no transitions in the histogram
xlabel(['almost T1 transitions'])

figure()
count_data = [hist(almostT1count,[1:4:30])]; %include the junctions with no transitions in the histogram
bar(count_data)
set(gca,'Xtick',[1:8],'XtickLabel',{'1','5','9','13','17','21','25','29'})
xlim([0.5 8.5])
ylim([0 20])
xlabel(['almost T1 transitions'])

%Histogram of the time before the same junction next contracts to a four
%way vertex
time_i = 0;
interval_i = 0;
almost_T1_interval = [];
stable_almost_T1s = [];
for i = unique(almost_T1_ids)
    time_i = almost_T1_time(find(almost_T1_ids == i)); 
    if length(time_i) > 1
        interval_i = time_i(2:end)-time_i(1:end-1);
        almost_T1_interval = [almost_T1_interval interval_i];
    else
        stable_almost_T1s = [stable_almost_T1s i];
    end
end

figure()
hist(almost_T1_interval)
xlabel('interval between transitions')

%Histogram of the number of T1 transitions per junction
T1count = 0;

[T1count,un]=hist(T1_ids,unique(T1_ids));

figure()
count_data = [hist(T1count,[1:1:8])]; %include the junctions with no transitions in the histogram
bar(count_data)
set(gca,'Xtick',[1:8],'XtickLabel',{'1','2','3','4','5','6','7','8'})
xlim([0.5 8.5])
ylim([0 18])
xlabel(['T1 transitions'])

%Histogram of the time before the same junction undergoes another T1
%transition

time_i = 0;
interval_i = 0;
T1_interval = [];
stable_T1s = [];
for i = unique(T1_ids)
    time_i = T1_time(find(T1_ids == i)); 
    if length(time_i) > 1
        interval_i = time_i(2:end)-time_i(1:end-1);
        T1_interval = [T1_interval interval_i];
    else
        stable_T1s = [stable_T1s i];
    end
end

figure()
interval_data = [hist(T1_interval/2,[4:8:120]) 0 0 0 0 16];
bar(interval_data)
set(gca,'Xtick',[1:20],'XtickLabel',{'4','12','20','28','36','42','','','68','','84','','','108','116','','','','','stable'})
xlim([0.5 20.5])
ylim([0 18])
xlabel('Interval between T1 events [minutes]')