big_table = importfile('Data.xlsx');    % Importing table with all the data

%%
c = 1;
close all
for ff = 1%:size(Events,2);
    volume = big_table(2:end,(ff-1)*3+1);    % Exctracting single cells volume
    GFP_conc = big_table(2:end,(ff-1)*3+2);  % Extracting singel cells GFP concentration
    
    set = volume > 0;   
    volume = volume(set);   % Elimintates NaNs value at the end
    GFP_conc = GFP_conc(set);   % Elimintates NaNs value at the end
    
    vol = volume;  %Extracting single cell cycles
    gfp_conc = GFP_conc;   %Extracting single cell cycles
    total = vol.*gfp_conc;  %calculate un-smoothed total GFP
    
    time = linspace(1, size(vol,1),size(vol,1))';
    
    % A smoothing spline will be used to smooth the total GFP
    
    %fitting total GFP with a smoothing spline
    f1 = fit(time,total,'smoothingspline','SmoothingParam',0.03); 
    smoothed_tot = f1(time);
    
    %differentiating the smoothed total GFP
    
    tot_diff = differentiate(f1,time);
    
    
    %fitting with a tight spline the first derivatives of the total GFP to differentiate it again
    
    tot_diff_smooth = fit(time,tot_diff,'smoothingspline','SmoothingParam',1); 
    tot_diff2 = differentiate(tot_diff_smooth,time);
    
    
    %calculate GFP production rate assuming a maturatio time(Tm) to be 6 min
    
    Km = 0.346;  %Km = ln2/Tm
    
    % calculate GFP production rate corrected for the GFP maturation time
    prod_rate = tot_diff + (1/Km)*tot_diff2;
    
    % calculate uncorrected GFP production rate
    prod_rate_u = tot_diff;
    
    % plot raw total GFP and and smoothing to visually check it
%     figure('pos',[50 50 1100 700])
%     subplot(1,2,1)
%     plot(time,total,'.','MarkerSize',12);
%     hold on
%     plot(f1);
%     xlim([time(1) time(end)])
%     subplot(1,2,2)
%     plot(time,prod_rate,'-b','MarkerSize',12);
%     xlim([time(1) time(end)])
       
    c = c+1;
  
end
%%

% save raw and smoothed variables
save(strcat('cell_cycle_',num2str(ff)),'smoothed_tot','tot_diff','tot_diff2','prod_rate','prod_rate_u','gfp_conc','vol','total');