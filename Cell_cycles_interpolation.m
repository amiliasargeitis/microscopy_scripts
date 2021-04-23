
list = dir('**.mat');   %create alist with the names of all the .mat files containing singl cell cycle variables
list2 = natsortfiles({list.name})'; %order the list of .mat files
figure,
%%
for n = 1:size(list2)
    
    table = load(list2{n}); %load single cell cycle variables
    time = linspace(1,size(table.vol,1),size(table.vol,1)); %define a real time axis, in frames
    
    bud = Events(1,n);  %extract budding event for the loaded cell cycle (starting from the begginning i.e. 3 frames before the first karyokinesis)
    
    time1 = linspace(1,bud-4,bud-4);    %define a time axis (1 to budding) for the first part of the cell cycle (from karyokinesis(4th frame) to the budding)
    time2 = linspace(bud-3,time(1,end)-7,time(1,end)-bud-3);    %define a time axis (budding to karyokinesis) for the second part of the cell cycle (from budding to the second karyokinesis((end-3)th frame))
    
    % define relative time axis for interpolation of the first part(G1) of the
    % cell cycle with 22 points (this number is decided based on the ratio between the avrage G1 duration and the average S/G2/M duration for the analyzed cell cycles)
    Time1_1 = time1;
    Time1_2 = [0:1/(size(Time1_1,2)-1):1];
    Time1_3 = [0:1/21:1];
    
    % define relative time axis for interpolation of the second part(S/G2/M) of the
    % cell cycle with 38 points (this number is decided based on the ratio between the avrage G1 duration and the average S/G2/M duration for the analyzed cell cycles)
    Time2_1 = time2;
    Time2_2 = [0:1/(size(Time2_1,2)-1):1];
    Time2_3 = [0:1/37:1];
    
    % interpolate each variable splitted in the first and second part of the cell cycle with a fixed
    % number of points
    rate1_1 = table.prod_rate(5:bud);
    rate1_2 = interp1(Time1_2,rate1_1',Time1_3,'spline');
    
    rate2_1 = table.prod_rate(bud+1:time(1,end)-3);
    rate2_2 = interp1(Time2_2,rate2_1',Time2_3,'spline');
    
    rate1_1u = table.prod_rate_u(5:bud);
    rate1_2u = interp1(Time1_2,rate1_1u',Time1_3,'spline');
    
    rate2_1u = table.prod_rate_u(bud+1:time(1,end)-3);
    rate2_2u = interp1(Time2_2,rate2_1u',Time2_3,'spline');
   
    % combine the two interpolated part for each variable of interest and
    % save them in a matrix
    RATES(:,n) = [rate1_2,rate2_2];
    RATESu(:,n) = [rate1_2u,rate2_2u];
    
end
