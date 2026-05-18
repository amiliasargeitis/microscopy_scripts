#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
Modules required for imaging and data analysis
'''

import pandas as pd
import numpy as np
import matplotlib##
import matplotlib.pyplot as plt##
import scipy.interpolate
get_ipython().run_line_magic('matplotlib', 'inline')
import os
import math
from math import pi, sin, cos
import seaborn as sns
import scipy
import skimage
import random
import re
from scipy import interpolate
from scipy.ndimage import binary_erosion
from skimage.io import imread, imshow
from scipy.ndimage import center_of_mass
from skimage import exposure, data
from skimage import io
import os
from skimage import filters
import math
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
from skimage.filters import threshold_local, threshold_otsu, threshold_multiotsu, threshold_li, threshold_yen, threshold_minimum
from skimage.morphology import remove_small_objects, disk, erosion
import statsmodels.api as sm
import bottleneck as bn
nanmean = bn.nanmean
import cv2 as cv
import sklearn


# In[14]:


'''
Conversion from TimeIDs to time points in minutes.

217 time points are in the movie.

'''
#### numbers different for different experiments:
total_number_of_TimeIDs = 109

#time_loop_before_the_switch = 95
#time_loop_after_the_switch = time_loop_before_the_switch + 1

time_step = 5

TimeIDs = range(1, total_number_of_TimeIDs + 1, 1) #+1 because the last number is not accounted
Time_mins = []

for i in TimeIDs:
    t = time_step*(i-1)
    Time_mins.append(t)
    
    
time_conversion = pd.DataFrame({"TimeID" : TimeIDs, "Time" : Time_mins})

time_conversion


# In[3]:


'''
!!!Not necessary to run each time!!!
'''

PATH_to_the_TIFs = "C:\\Users\\andri\\Documents\\PhD\\Temporary folder\\RF2_SFP1_Light\\Processed_tiffs\\"
all_files = os.listdir(PATH_to_the_TIFs)

tif_files = []
for f in all_files:
    if f.endswith(".tif"):
        tif_files.append(f)
        
tif_files


# In[4]:


max_time = list(time_conversion["Time"])[-1]
min_time = list(time_conversion["Time"])[0]

print (max_time, min_time)


# In[5]:



'''
Creating a list of files from BudJ to open plus prefixes
'''
PATH_to_the_TIFs2 = "C:\\Users\\andri\\Documents\\PhD\\Temporary folder\\RF2_SFP1_Light\\Analysis\\Budj_tracking\\"
files = []
for i in os.listdir(PATH_to_the_TIFs2):
    if ".csv" in i:
            prefix_position = "pos" + i[-25:-23]+"_"
            files.append((i, prefix_position))
            
files


# In[6]:


'''
Creating a DataFrame containing the segmenation info from BudJ
'''

initial_table = pd.DataFrame({})
for f in files:
    pos = pd.read_csv(PATH_to_the_TIFs2+f[0], header=0, index_col=0)
    pos["Cell_pos"] = f[1] + pos["Cell"].map(str)
    f_columns = ["TimeID", "Cell_pos", "Volume", "x", "y", "Major R", "Minor r", "Angle"]
    pos = pos.reindex(columns = f_columns)
    initial_table = pd.concat([initial_table, pos])
    
initial_table = pd.merge(initial_table, time_conversion, on="TimeID")
initial_table


# In[7]:


'''
Creating a a list of files from Click_cells 
'''

PATH_to_the_macro_clicking = "C:\\Users\\andri\\Documents\\PhD\\Temporary folder\\RF2_SFP1_Light\\Analysis\\Events_daughters\\"
files2 = []
for i in os.listdir(PATH_to_the_macro_clicking):
    if ".csv" in i:
            prefix_position = "pos" + i[-6:-4]+"_"
            files2.append((i, prefix_position))
files2


# In[29]:


time_conversion = pd.DataFrame({"t" : TimeIDs, "Time" : Time_mins})

max_time = list(time_conversion["Time"])[-1]
min_time = list(time_conversion["Time"])[0]

initial_table2 = pd.DataFrame({})

for f in files2:
    
    pos1 = pd.read_csv(PATH_to_the_macro_clicking+f[0], header=0, index_col=0)
    pos1["Cell_pos"] = f[1] + pos1["Cell"].map(str)
    d_columns = ["Event", "Cell_pos", "x", "y", "t"]
    pos1 = pos1.reindex(columns = d_columns)
    initial_table2 = pd.concat([initial_table2, pos1])
    
new = pd.merge(initial_table2, time_conversion, on="t").sort_values("Time")
new


# In[28]:


'''
Creates a dictionary containing the birth event (in frames) of each cell 
'''

birth_Sfp1 = new.loc[new['Event'] == 'Start'].groupby('Cell_pos')['t'].apply(list).to_dict()


# In[37]:


'''
Creates a dictionary containing the budding events (in frames) of each cell 
'''

buddings_Sfp1 = (new.loc[new['Event'] == 'Bud'].groupby('Cell_pos')['t'].apply(list).to_dict())


# In[38]:


'''
Creates a dictionary containing the cytokinesis events (in frames) of each cell 
'''

cytokinesis_Sfp1 = (new.loc[new['Event'] == 'Cytokinesis'].groupby('Cell_pos')['t'].apply(list).to_dict())


# In[18]:



individual_cells2 = sorted(list(set(initial_table["Cell_pos"])))
print(len(individual_cells2)
     )


# In[19]:


'''
!!!Not necessary to run each time!!!
'''

#function to extract parameters of the ellipse from the BudJ table 
def ellipse(time_point, BudJ_table, scaling_factor):
    h = float(BudJ_table[BudJ_table["TimeID"] == time_point]['x'])/scaling_factor
    k = float(BudJ_table[BudJ_table["TimeID"] == time_point]['y'])/scaling_factor
    a = float(BudJ_table[BudJ_table["TimeID"] == time_point]["Major R"])/scaling_factor
    b = float(BudJ_table[BudJ_table["TimeID"] == time_point]["Minor r"])/scaling_factor
    A = float(BudJ_table[BudJ_table["TimeID"] == time_point]['Angle'])*(math.pi/180)
    return h, k, a, b, A


# In[20]:


def round_up_to_odd(f):
    return int(np.ceil(f) // 2 * 2 + 1)


# In[64]:


scaling_factor = 0.16 #microns per pixel, 100x objective

bloc_size_frac_to_use = 0.05
#offset_to_use = -15


initial_table_recalculated = pd.DataFrame({})

c = 0
for pos in ['05','06', '08', '09', '11','16','18','21', '22', '23', '24', '25', '26', '27', '28', '29', '31', '32', '33', '34', '37', '39', '40']:

        filename = os.path.join(skimage.data_dir, PATH_to_the_TIFs+"rf2+sfp1-mscarlet_light_xy"+pos+".nd2.tif")
        print(filename)
        
        cell_tif = io.imread(filename)

    
        for cell in individual_cells2:
            if 'pos'+pos in cell:
                nothing_found = 0
                temp = initial_table[initial_table["Cell_pos"] == cell]
                temp = temp.sort_values(by="Time")
                thresholded_frames = []
                time_axis = list(temp["TimeID"])
                
                
                #GFP_total = []
                RFP_cyto = []
                RFP_nucleus = []
                Ratio = []
                RFP_total = []
                RFP_std = []
                
                              
                for t in time_axis:
                    
                        
                        t_in_tiff = int(t) - 1

                        h, k, a, b, A = ellipse(t, temp, scaling_factor)

                        #creating the mask corresponding to the BudJ ellipse
                        image = cell_tif[t_in_tiff,:,:,2] #does not matter which channel
                        nrows, ncols = image.shape
                        row, col = np.ogrid[:nrows, :ncols]

                        inner_disk_mask = ((((col-h-1)*math.cos(A)+(row-k-1)*math.sin(A))**2)/(a**2) + (((col-h-1)*math.sin(A)-(row-k-1)*math.cos(A))**2)/(b**2) - 1 < 0)

       
                        #adaptive thresholding for the nucleus
                        imageRFP = cell_tif[t_in_tiff,:,:,2]
                        imageRFP_1 = cv.GaussianBlur(imageRFP, (3,3), cv.BORDER_DEFAULT)
                        nrows, ncols = imageRFP_1.shape
                        
                        row, col = np.ogrid[:nrows, :ncols]

                        inner_disk_mask = ((((col-h)*math.cos(A)+(row-k)*math.sin(A))**2)/(a**2) + (((col-h)*math.sin(A)-(row-k)*math.cos(A))**2)/(b**2) - 1 < 0)

                        image_masked_flattened = imageRFP_1[inner_disk_mask == True]

                        num_cell_pixels = len(image_masked_flattened)
                        #print(num_cell_pixels)
                        
                        bloc_size_cell_size_dependent = round_up_to_odd(bloc_size_frac_to_use*num_cell_pixels)


                        #adaptive thresholding will be performed not on all image but on a fraction of it which includes the cell
                        imageRFP_corr = imageRFP_1*inner_disk_mask
                    
                        #Adjust local threshold offset until nucleus area ≤ 30% of cell area.
                        
                        num_nuc_pixels = num_cell_pixels
    
                        offset_to_use = 5
                        while num_nuc_pixels > 0.3*num_cell_pixels: #nucleus corresponds to ~30% of the total cell area 
                            nuc_thresh_mask_local = threshold_local(image=imageRFP_corr, block_size=bloc_size_cell_size_dependent, offset=offset_to_use)
                            nuc_thresh_mask_local = imageRFP_corr > nuc_thresh_mask_local
                         
                            nuc_thresh_mask_local_flattened = imageRFP_1[nuc_thresh_mask_local == True]
                            num_nuc_pixels = len(nuc_thresh_mask_local_flattened)
                        
                            offset_to_use = offset_to_use-3 
                            continue
                            
                        # len(nuc_thresh_mask_local)
                        nuc_thresh_mask_local = remove_small_objects(nuc_thresh_mask_local, min_size=10)
                        nuc_thresh_mask_local_eroded = binary_erosion(nuc_thresh_mask_local, structure = np.ones((3,3)))
                        
                        
                        centroid = scipy.ndimage.center_of_mass(nuc_thresh_mask_local_eroded)
                        centroid = np.nan_to_num(centroid)
                        b = centroid[1]
                        a = centroid[0]
                        r = 3
                        r1 = 11
                            
                        x1,y1 = np.ogrid[-a: nrows-a, -b: ncols-b]
                        
                        disk_mask_nuc = x1*x1 +y1*y1 < r*r
                        disk_mask_cyto = x1*x1 +y1*y1 < r1*r1
                    
                    
                        #RFP
                        mean_RFP = np.mean(image_masked_flattened)
                        if mean_RFP < 5:
                            mean_RFP = np.nan
                        RFP_total.append(mean_RFP)
                        RFP_std.append(np.std(imageRFP[inner_disk_mask == True]))
                    
                        #RFP in the nucleus
                        nucleus_mean = np.mean(imageRFP_corr[disk_mask_nuc == True])
                        if nucleus_mean < 5:
                            nucleus_mean = np.nan
                        RFP_nucleus.append(nucleus_mean)
                    
                        #RFP in the cytosol
                        diff = np.logical_and(disk_mask_cyto, inner_disk_mask)
                        mask_of_cytoplasm = inner_disk_mask ^ diff
                        
                        
                        image_cyto = imageRFP*mask_of_cytoplasm
                        cyto_mean = np.mean(imageRFP[mask_of_cytoplasm == True]) 
                        if cyto_mean < 5:
                            cyto_mean = np.nan
                        RFP_cyto.append(cyto_mean)
                    
                        #Nuclear to cytosolic ratio
                        Ratio.append((nucleus_mean)/(cyto_mean))
                    
               # if there's no thresholded pixels, then show the masked whole cell (makes it easy to recognize where it went 'well' and where it did not)
                        if np.count_nonzero(nuc_thresh_mask_local_eroded) == 0:
                            nothing_found += 1
                        thresholded_frames.append((t, h, k, imageRFP * inner_disk_mask, a, b, imageRFP_1*(mask_of_cytoplasm +disk_mask_nuc)))
                print(f"{cell} - No threshold result {round(nothing_found / len(temp['TimeID']) * 100)} % of the time")

                temp["RFP_nucleus"] = pd.Series(RFP_nucleus, index=temp.index)
                temp["RFP_cyto"] = pd.Series(RFP_cyto, index=temp.index)
                temp["RFP_total"] = pd.Series(RFP_total, index=temp.index)
                temp["RFP_std"] = pd.Series(RFP_std, index=temp.index)
                temp["Ratio"] = pd.Series(Ratio, index=temp.index)
                
             
                initial_table_recalculated = pd.concat([initial_table_recalculated, temp])
                
                c += 1
                print (cell,)


# In[65]:


initial_table_recalculated = initial_table_recalculated.reset_index(drop="true")


# In[66]:


initial_table_recalculated[pd.isnull(initial_table_recalculated).any(axis=1)]


# In[22]:


initial_table_recalculated.to_excel("Sfp1_mScarlet_RF2_light_2.xlsx")


# In[21]:


initial_table_recalculated = pd.read_excel("Sfp1_mScarlet_RF2_light_2.xlsx")


# In[22]:


initial_table = initial_table_recalculated
initial_table


# In[23]:


''''Extracting single cells'''

start_matrix = pd.DataFrame({})
buddings_matrix = pd.DataFrame({})
cytokinesis_matrix = pd.DataFrame({})

for cell in individual_cells2:

    start_matrix_single =  pd.DataFrame({})
    buddings_matrix_single = pd.DataFrame({})
    cytokinesis_matrix_single = pd.DataFrame({})
    table = initial_table[initial_table.Cell_pos == cell]
    table.to_excel(cell+".xlsx")
    start_matrix_single[cell] = birth_Sfp1[cell]
    buddings_matrix_single[cell] = buddings_Sfp1[cell]
    cytokinesis_matrix_single[cell] = cytokinesis_Sfp1[cell]
    start_matrix = pd.concat([start_matrix, start_matrix_single], axis = 1)
    buddings_matrix = pd.concat([buddings_matrix,buddings_matrix_single], axis = 1)
    cytokinesis_matrix = pd.concat([cytokinesis_matrix,cytokinesis_matrix_single], axis = 1)
    
start_matrix.to_csv("Start.csv")    
buddings_matrix.to_csv("Buddings.csv")  
cytokinesis_matrix.to_csv("Cytokinesis.csv")  


# In[39]:


start_control_SFP1 = {}


for cell in individual_cells2:
        control = []
        for i in birth_Sfp1[cell]:
            start_time = float(time_conversion[time_conversion["t"] == i]["Time"])
            control.append(i)
            
            
        start_control_SFP1[cell] = control


# In[27]:


buddings_control_SFP1 = {}


for cell in individual_cells2:
        control = []
        for i in buddings_Sfp1_daughters[cell]:
            budding_time = float(time_conversion[time_conversion["t"] == i]["Time"])
            control.append(i)
            print(control)
            
        buddings_control_SFP1[cell] = control


# In[28]:


cytokinesis_control_SFP1 = {}

 
for cell in individual_cells2:
        control = []
        for i in cytokinesis_Sfp1_daughters[cell]:
            cytokinesis_time = float(time_conversion[time_conversion["t"] == i]["Time"])
            control.append(i)
            print(control)
            print(i)
        cytokinesis_control_SFP1[cell] = control
        print(cytokinesis_control_SFP1)


# In[29]:


'''alignment cytokinesis to cytokinesis'''
'''for daughter cells using Events: Start, Budding, Cytokinesis'''

cell_cycles_list_SFP1_d = []

start_group_SFP1 = start_control_SFP1
buddings_group_SFP1 = buddings_control_SFP1
cytokinesis_group_SFP1 = cytokinesis_control_SFP1
flag = "CONTROL"


for cell in cytokinesis_group_SFP1:
    
    table = initial_table_recalculated[initial_table_recalculated["Cell_pos"] == cell]
    
    birth_of_the_cell = start_group_SFP1[cell]
    buddings_of_the_cell = buddings_group_SFP1[cell]
    cytokinesis_of_the_cell = cytokinesis_group_SFP1[cell]
    print(len(cytokinesis_of_the_cell))
    
    for i in range(len(birth_of_the_cell)):
        #if i != len(cytokinesis_of_the_cell)-1:
            print (cell)
            
            start1 = float(time_conversion[time_conversion["t"] == birth_of_the_cell[i]]["Time"])
            #print(start1)
            start1 = start1 -10 
            end1 = float(time_conversion[time_conversion["t"] == buddings_of_the_cell[i]]["Time"])### KEEP the FLOAT function here
            #print(end1)
            start2 = float(time_conversion[time_conversion["t"] == buddings_of_the_cell[i]]["Time"])
            #print(start2)
            end2 = float(time_conversion[time_conversion["t"] == cytokinesis_of_the_cell[i]]["Time"])### KEEP the FLOAT function here
            #print(end2) 
            end2 = end2 -10
            
            duration1 = (end1 - start1) #G1 duration of daughter cells
            #print(duration1)
        
            duration2 = (end2 - start2) #SG2M duration of daughter cells 
            #print(duration2)
            time_points_for_aligning_at_bud1 = (1/duration1)*np.array(np.arange(start1-start1, end1-start1+1, 5))
            time_points_for_aligning_at_bud2 = (1/duration2)*np.array(np.arange(start2-start2, end2-start2+1, 5))
            #print(time_points_for_aligning_at_bud1)
            #print(time_points_for_aligning_at_bud2)
            cell_cycles_list_SFP1_d.append((cell, start1, end1, start2, end2, time_points_for_aligning_at_bud1, time_points_for_aligning_at_bud2, duration1, duration2))
            print(cell_cycles_list_SFP1_d)
 
            
            


# In[30]:


import matplotlib.pyplot as plt

# Data and plot setup (assuming 'cell_cycles_list_SFP1_d' is defined)
sg2m_d = []
g1_d = []
tot_d = []

for cycles in cell_cycles_list_SFP1_d:
    dur1 = cycles[7]
    dur2 = cycles[8]
    dur3 = cycles[7] + cycles[8]
    sg2m_d.append(dur2)
    g1_d.append(dur1)
    tot_d.append(dur3)

# Plotting
plt.figure(figsize=(9, 6))

plt.hist(tot_d, bins=15, alpha=0.5, label='Full cell cycle duration', color='purple')
plt.hist(g1_d, bins=15, alpha=0.5, label='G1 duration', color='orange')
plt.hist(sg2m_d, bins=15, alpha=0.5, label='S/G2/M duration', color='green')

plt.xlabel("Duration (min)", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.title("Distribution of Cell Cycle Durations (Daughters)", fontsize=14)
plt.legend(fontsize=10)
plt.xlim(0, 300)

plt.savefig("cell_cycle_hist_RF2_daughters_light.png")
plt.show()


# In[31]:


print(np.mean(tot_d))
print(np.mean(g1_d))
print(np.std(g1_d))
print(np.mean(sg2m_d))
print(np.mean(g1_d)/np.mean(tot_d))


# In[34]:


variables = ["Volume","RFP_cyto","RFP_nucleus","Ratio","RFP_total", "RFP_std"]

L = len(variables)
k=0

c = 1

all_time_points1 = np.linspace(0, 1, 80) #total number of points used for the interpolation of each cell cycle series
all_time_points2 = np.linspace(0, 1, 31) #numnber of points used to interpolate the first part (G1) of each cell cycle series
all_time_points3 = np.linspace(0, 1, 50) #numnber of points used to interpolate the second part (SG2M) of each cell cycle series
                                         #+1 point that will be exclude since it is already contained in the fisrt part (G1) of the cell cycle



cc_df1_SFP1_d = pd.DataFrame({})
cc_df_small1_SFP1_d = pd.DataFrame({})

cell_cycle = 0

for cycle in cell_cycles_list_SFP1_d:
    pos_cell = cycle[0]
    start1, end1, start2, end2, time1, time2 = cycle[1:7]
    new_time1, new_time2, new_time3 = all_time_points1, all_time_points2, all_time_points3
    cell_cycle += 1
    print(f"\n--- Processing Cell {pos_cell} ---")

    # reset each cycle
    cc_df_small1_SFP1_d = pd.DataFrame({})

    for variable in variables:
        big_data_table = initial_table_recalculated
        table = big_data_table[big_data_table["Cell_pos"] == pos_cell]

        sensor1 = table[(table["Time"] >= start1) & (table["Time"] <= end1)][variable]
        sensor2 = table[(table["Time"] >= start2) & (table["Time"] <= end2)][variable]

        # check for missing values
        if sensor1.empty or sensor2.empty:
            print(f" {variable}: No data in selected time ranges "
                  f"(len(sensor1)={len(sensor1)}, len(sensor2)={len(sensor2)})")
            continue

        # check for length mismatch
        if len(sensor1) != len(time1):
            print(f"{variable}: len(sensor1)={len(sensor1)} != len(time1)={len(time1)}")
            continue
        if len(sensor2) != len(time2):
            print(f"{variable}: len(sensor2)={len(sensor2)} != len(time2)={len(time2)}")
            continue

        try:
            # linear interpolation
            f1 = scipy.interpolate.interp1d(time1, sensor1, bounds_error=False, fill_value=np.nan)
            f2 = scipy.interpolate.interp1d(time2, sensor2, bounds_error=False, fill_value=np.nan)

            new_sensor1 = f1(new_time2)
            new_sensor2 = f2(new_time3)

            
            if np.all(np.isnan(new_sensor1)) or np.all(np.isnan(new_sensor2)):
                print(f"{variable}: interpolation gave only NaNs")
                continue

           
            cc_df_small1_SFP1_d["TimeID"] = new_time1
            cc_df_small1_SFP1_d[variable] = np.append(new_sensor1, new_sensor2[1:])
            print(f"{variable}: interpolated successfully "
                  f"(len={len(new_sensor1)+len(new_sensor2)-1})")

        except Exception as e:
            print(f" {variable}: interpolation failed with error: {e}")
            continue

    # only add if something was filled
    if not cc_df_small1_SFP1_d.empty:
        cc_df_small1_SFP1_d["Cell_cycle"] = pos_cell
        cc_df_small1_SFP1_d["Ratio3"] = cell_cycle
        cc_df1_SFP1_d = pd.concat([cc_df1_SFP1_d, cc_df_small1_SFP1_d], ignore_index=True)
    else:
        print(f"Skipped cell {pos_cell}: no valid variables added")
                


# In[36]:


cc_df1_SFP1_d.to_excel("20240130_RF2_light_daughters2.xlsx")
cc_df1_SFP1_d = pd.read_excel("20240130_RF2_light_daughters2.xlsx")


# In[37]:


'''
Normalize each single-cell trace using Min-Max scaling to the range [-1, 1]
'''

from sklearn.preprocessing import MinMaxScaler

scaler=MinMaxScaler(feature_range=(-1, 1))

orig = cc_df1_SFP1_d.pivot(index='Ratio3', columns='TimeID', values='Ratio')

scaled=pd.DataFrame(scaler.fit_transform(orig.T).T,columns=orig.columns)


# In[40]:


scaled.to_csv("OC_3_LIGHT_heatmap.csv")


# In[41]:


variables = ["Volume","RFP_total","RFP_std", "RFP_nucleus","RFP_cyto","Ratio"]
ynames = ["Volume","RFP_total","RFP_std","RFP_nucleus","RFP_cyto","Ratio"]
colors = ["grey","red","red","red","red","blue"]


L = len(variables)

plt.ioff()
fig = plt.figure(1, (5*L, 5))
sns.set_style("darkgrid")
c = 1

all_time_points = np.linspace(0, 1, 80)
        
for i in range(L):
    variable = variables[i]
    ax = plt.subplot(1, L, c)
    sns.tsplot(data=cc_df1_SFP1_d, time="TimeID", unit="Cell_cycle", value=variable, color=colors[i], estimator=nanmean)
    
    c+=1
    plt.xlabel("Cell life cycle")
    plt.ylabel(ynames[i])
    #plt.ylim(0.8,1.1)
    plt.title("82 cells,82 cycles ")
    ax.axvline(0.4, color='black', linewidth=1, linestyle='--')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(14)
    
plt.tight_layout()
plt.savefig("cc_aligned_RF2_daughters_light2.png")
#plt.savefig("comparison_SFP1_TOD6_b2b.png")
plt.show()

