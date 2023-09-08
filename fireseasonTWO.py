import numpy as np
import netCDF4 as nc
from datetime import datetime

####################################################################################################
# LOAD IN NETCDF FILES
# temp, humidity, precip
fn1 = 'C:/Users/modarn/Documents/CASQA/data/ta_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc'
fn2 = 'C:/Users/modarn/Documents/CASQA/data/hur_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc'
fn3 = 'C:/Users/modarn/Documents/CASQA/data/pr_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc'
data1IN = nc.Dataset(fn1)
data2IN = nc.Dataset(fn2)
data3IN = nc.Dataset(fn3)
## define varaibles from netCDF files
# tempK = data1IN['ta'][:]
tempC_in = data1IN['ta'][:] - 272.16
RH_in = data2IN['hur'][:]
prec_in = data3IN['pr'][:]
prec_in_inches = prec_in * 3401.2 #in/day
# 1kg = 1Liter = 10e6 mm3 and 1m2 = 10e6 mm2. 1kg of rain over 1 m2 is equivalent to 1mm
# 1 kg/m2/s = 86400 mm/day  factor of 86400 (seconds in a day)
lat = data1IN['lat'][:]
lon = data1IN['lon'][:]
## TIME - SEASON 2
# 1997/1998 season, fire season = May through August
# Extract the time data and convert it to datetime objects
time_data = data1IN.variables['time'][:]
time_units = data1IN.variables['time'].units
time_calendar = data1IN.variables['time'].calendar
time_datetime = nc.num2date(time_data, units=time_units, calendar=time_calendar)
# Fire season - define time
# Define the date range from May 1, 1997, to August 28, 1997 (sept. 1st but reducing days for data compatibility)
start_date = datetime(1998, 5, 1)
end_date = datetime(1998, 9, 28)
# Find the indices of time data that fall within the specified date range
time_indices = [idx for idx, t in enumerate(time_datetime) if start_date <= t <= end_date]
# Create a variable that shows the datetime values selected by time_indices to check work
selected_datetime_values = [time_datetime[idx] for idx in time_indices]
# Wet season - define time
# Define the date range from November 1, 1997 to March 31, 1998
start_date_wet = datetime(1997, 11, 1)
end_date_wet = datetime(1998, 3, 31)
# Find the indices of time data that fall within the specified date range
time_indices_wet = [idx for idx, t in enumerate(time_datetime) if start_date_wet <= t <= end_date_wet]
# Create a variable that shows the datetime values selected by time_indices to check work
selected_datetime_values_wet = [time_datetime[idx] for idx in time_indices_wet]
# Extract the data for the specified date range for PRECIP, TEMPERATURE AND RELATIVE HUMIDITY VARIABLE
tempC = tempC_in[time_indices, :, :, :]  # SAVE THESE VALUES
RH = RH_in[time_indices, :, :, :]
prec = prec_in_inches[time_indices_wet, :, :]
# CALCULATE DEWPOINT TEMPERATURE
# Initialize empty matrices to store the dew point temperatures for low and high elevations
Td_low = np.zeros_like(tempC[:, 1, :, :])
Td_high = np.zeros_like(tempC[:, 2, :, :])
# Calculate the total number of iterations for each elevation
total_iterations_low = tempC.shape[0] * tempC.shape[2] * tempC.shape[3]
total_iterations_high = tempC.shape[0] * tempC.shape[2] * tempC.shape[3]
# Initialize counters for the iterations
iteration_count_low = 0
iteration_count_high = 0
# Loop through each time, latitude, and longitude to calculate the dew point temperature
for t in range(tempC.shape[0]):
    for lat in range(tempC.shape[2]):
        for lon in range(tempC.shape[3]):
            # Calculate dew point temperature using the provided formula for low elevation
            T_low = tempC[t, 1, lat, lon]
            rh_low = RH[t, 1, lat, lon] / 100.0
            Td_low[t, lat, lon] = (237.3 * (np.log(rh_low) + (17.27 * T_low) / (237.3 + T_low))) / \
                                  (17.27 - (np.log(rh_low) + (17.27 * T_low) / (237.3 + T_low)))

            # Increment the low elevation iteration count and print progress
            iteration_count_low += 500000
            # Print the dew point temperature for each time, latitude, and longitude
            # print(f"Dew Point Temperature at time {t}, lat {lat}, lon {lon}: {Td_low[t, lat, lon]}")
            # print(f"Low Elevation - Iteration {iteration_count_low}/{total_iterations_low}")

            # Calculate dew point temperature using the provided formula for high elevation
            T_high = tempC[t, 2, lat, lon]
            rh_high = RH[t, 2, lat, lon] / 100.0
            Td_high[t, lat, lon] = (237.3 * (np.log(rh_high) + (17.27 * T_high) / (237.3 + T_high))) / \
                                   (17.27 - (np.log(rh_high) + (17.27 * T_high) / (237.3 + T_high)))

            # Increment the high elevation iteration count and print progress
            iteration_count_high += 1
            # print(f"Dew Point Temperature at time {t}, lat {lat}, lon {lon}: {Td_high[t, lat, lon]}")
            # print(f"High Elevation - Iteration {iteration_count_high}/{total_iterations_high}")
# Now the dew point temperatures for low and high elevations are stored in 'Td_low' and 'Td_high' respectively
########################################################################################
# SEPERATE DATA BASED ON ELEVATION LEVELS (TANG ET AL., 2015)
# LAPSE RATE (component(A) - Stability component based on elevation variants
lowA = tempC[:, 0, :, :] - tempC[:, 1, :, :]
midA = tempC[:, 1, :, :] - tempC[:, 2, :, :]
highA = tempC[:, 2, :, :] - tempC[:, 3, :, :]
#####
# DEW POINT (component(B)) - humidity component based on elevation variants
lowB = tempC[:, 1, :, :] - Td_low[:, :, :]
midB = tempC[:, 1, :, :] - Td_low[:, :, :]
highB = tempC[:, 2, :, :] - Td_high[:, :, :]
#########################################################################
# assign HI index to values
# Create a new array to store the assigned numeric constants
lowA_HI = np.zeros_like(lowA)
midA_HI = np.zeros_like(midA)
highA_HI = np.zeros_like(highA)
lowB_HI = np.zeros_like(lowB)
midB_HI = np.zeros_like(midB)
highB_HI = np.zeros_like(highB)
## Assign numeric constants based on temperature values
# Low A
lowA_HI[np.where(lowA < 4)] = 1
lowA_HI[np.where((lowA >= 4) & (lowA <= 7))] = 2
lowA_HI[np.where(lowA >= 8)] = 3
# mid A
midA_HI[np.where(midA < 6)] = 1
midA_HI[np.where((midA >= 6) & (midA <= 10))] = 2
midA_HI[np.where(midA >= 11)] = 3
# high A
highA_HI[np.where(highA < 18)] = 1
highA_HI[np.where((highA >= 18) & (highA <= 21))] = 2
highA_HI[np.where(highA >= 22)] = 3
##########################################
# Low B
lowB_HI[np.where(lowB < 6)] = 1
lowB_HI[np.where((lowB >= 6) & (lowA <= 9))] = 2
lowB_HI[np.where(lowB >= 10)] = 3
# mid A
midB_HI[np.where(midB < 6)] = 1
midA_HI[np.where((midB >= 6) & (midA <= 12))] = 2
midA_HI[np.where(midB >= 13)] = 3
# high A
highB_HI[np.where(highB < 15)] = 1
highB_HI[np.where((highB >= 15) & (highA <= 20))] = 2
highB_HI[np.where(highB >= 21)] = 3
#
##################################################################
# FINAL HI VALUES
HI_low = (lowA_HI + lowB_HI)  # SAVE THIS VARIABLE
HI_mid = (midA_HI + lowB_HI)  # SAVE THIS VARIABLE
HI_high = (highA_HI + highB_HI)  # SAVE THIS VARIABLE
# Calculate the average across all time steps for HI_low
HI_low_average_time = np.mean(HI_low, axis=0)
HI_mid_average_time = np.mean(HI_mid, axis=0)
HI_high_average_time = np.mean(HI_high, axis=0)
# Round the average array up to the nearest integer
HI_low_avg = np.ceil(HI_low_average_time)  # SAVE THIS VARIABLE
HI_mid_avg = np.ceil(HI_mid_average_time)  # SAVE THIS VARIABLE
HI_high_avg = np.ceil(HI_high_average_time)  # SAVE THIS VARIABLE
## calculate the percentage of days greater than or equal to HI=5
# Count the number of days where HI values are greater than or equal to 5 for each category.
# Calculate the total number of days in each category.
# Divide the count of days greater than or equal to 5 by the total number of days in that category and multiply by 100 to get the percentage.
# Step 1: Count the number of days greater than or equal to HI = 5 for each category
HI_low_ge_5 = np.sum(HI_low >= 5)
HI_mid_ge_5 = np.sum(HI_mid >= 5)
HI_high_ge_5 = np.sum(HI_high >= 5)
# Step 2: Calculate the total number of days in each category
total_days_low = HI_low.size
total_days_mid = HI_mid.size
total_days_high = HI_high.size
# Step 3: Calculate the percentage of days greater than or equal to HI = 5 for each category
HI_low_perct = (HI_low_ge_5 / total_days_low) * 100
HI_mid_perct = (HI_mid_ge_5 / total_days_mid) * 100
HI_high_perct = (HI_high_ge_5 / total_days_high) * 100
# Print the results
# print(f"Percentage of days with HI greater than or equal to 5 in the low category: {HI_low_perct:.2f}%")
# print(f"Percentage of days with HI greater than or equal to 5 in the mid category: {HI_mid_perct:.2f}%")
# print(f"Percentage of days with HI greater than or equal to 5 in the high category: {HI_high_perct:.2f}%")
###########################################################################################################
# SPATIAL CORRELATION
# Calculate the spatial correlation between HI_low and prec for each grid point
spatial_correlations_HIlow = np.zeros_like(HI_low_avg)  # Initialize an array to store correlations
# Loop through each latitude and longitude
for lat in range(HI_low_avg.shape[0]):
    for lon in range(HI_low_avg.shape[1]):
        # Exclude points where either HI_low or prec has zero variance
        if np.var(HI_low[:, lat, lon]) != 0 and np.var(prec[:, lat, lon]) != 0:
            correlation_matrix = np.corrcoef(HI_low[:, lat, lon], prec[:, lat, lon])
            spatial_correlation_HIlow = correlation_matrix[0, 1]
            spatial_correlations_HIlow[lat, lon] = spatial_correlation_HIlow
# Calculate the spatial correlation between HI_mid and prec for each grid point
spatial_correlations_HImid = np.zeros_like(HI_mid_avg)  # Initialize an array to store correlations
# Loop through each latitude and longitude
for lat in range(HI_mid_avg.shape[0]):
    for lon in range(HI_mid_avg.shape[1]):
        # Exclude points where either HI_mid or prec has zero variance
        if np.var(HI_mid[:, lat, lon]) != 0 and np.var(prec[:, lat, lon]) != 0:
            correlation_matrix = np.corrcoef(HI_mid[:, lat, lon], prec[:, lat, lon])
            spatial_correlation_HImid = correlation_matrix[0, 1]
            spatial_correlations_HImid[lat, lon] = spatial_correlation_HImid
# Calculate the spatial correlation between HI_mid and prec for each grid point
spatial_correlations_HIhigh = np.zeros_like(HI_mid_avg)  # Initialize an array to store correlations
# Loop through each latitude and longitude
for lat in range(HI_high_avg.shape[0]):
    for lon in range(HI_high_avg.shape[1]):
        # Exclude points where either HI_high or prec has zero variance
        if np.var(HI_high[:, lat, lon]) != 0 and np.var(prec[:, lat, lon]) != 0:
            correlation_matrix = np.corrcoef(HI_high[:, lat, lon], prec[:, lat, lon])
            spatial_correlation_HIhigh = correlation_matrix[0, 1]
            spatial_correlations_HIhigh[lat, lon] = spatial_correlation_HIhigh
#######################################################################################################
# Create a NetCDF file to save the 'prec' variable
output_file_path = 'C:/Users/modarn/Documents/CASQA/data/season1997_1998.nc'
output_ncfile = nc.Dataset(output_file_path, 'w')
##### DEFINE EACH VARIABLE DIMENSIONS
#####################
# PRECIPITATION
# Define dimensions based on the shape of 'prec' variable
time_dim = output_ncfile.createDimension('time_wet', prec.shape[0])
lat_dim = output_ncfile.createDimension('latitude', prec.shape[1])
lon_dim = output_ncfile.createDimension('longitude', prec.shape[2])
# Create a NetCDF variable for 'prec' using appropriate dimensions and data type
nc_prec = output_ncfile.createVariable('precipitation', np.float32, ('time_wet', 'latitude', 'longitude'))
# Assign data from 'prec' to the NetCDF variable
nc_prec[:] = prec
# Add units and description attributes to the NetCDF variable if needed
nc_prec.units = 'inches'
nc_prec.description = 'Daily precipitation in inches'
del time_dim, lat_dim, lon_dim
######################
# TEMPERATURE
# Define dimensions based on the shape of 'tempC' variable.
time_dim = output_ncfile.createDimension('time', tempC.shape[0])
plev_dim = output_ncfile.createDimension('pressure_level', tempC.shape[1])
# Create a NetCDF variable for 'temp' using appropriate dimensions and data type
nc_tempC = output_ncfile.createVariable('temperature', np.float32, ('time','pressure_level', 'latitude', 'longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_tempC[:] = tempC
# Add units and description attributes to the NetCDF variable if needed
nc_tempC.units = 'Celcius'
nc_tempC.description = 'Daily avg temperature in Celcius'
######################
# RELATIVE HUMIDITY
# Create a NetCDF variable for 'RH' using appropriate dimensions and data type
nc_RH = output_ncfile.createVariable('relative_humidity', np.float32, ('time','pressure_level', 'latitude', 'longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_RH[:] = RH
# Add units and description attributes to the NetCDF variable if needed
nc_RH.units = 'percent'
nc_RH.description = 'Daily avg relative humidity values'
######################
# HI_low
nc_HIlow = output_ncfile.createVariable('HI_low',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HIlow[:] = HI_low
# Add units and description attributes to the NetCDF variable if needed
nc_HIlow.units = 'haines index'
nc_HIlow.description = 'HI values - low pressure level'
######################
# HI_mid
nc_HImid = output_ncfile.createVariable('HI_mid',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HImid[:] = HI_mid
# Add units and description attributes to the NetCDF variable if needed
nc_HImid.units = 'haines index'
nc_HImid.description = 'HI values - mid pressure level'
######################
# HI_high
nc_HIhigh = output_ncfile.createVariable('HI_high',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HIhigh[:] = HI_high
# Add units and description attributes to the NetCDF variable if needed
nc_HIhigh.units = 'haines index'
nc_HIhigh.description = 'HI values - high pressure level'
######################
# HI_avg_low
nc_HI_low_avg = output_ncfile.createVariable('HI_low_avg',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HI_low_avg[:] = HI_low_avg
# Add units and description attributes to the NetCDF variable if needed
nc_HI_low_avg.units = 'haines index'
nc_HI_low_avg.description = 'avg HI values - low pressure level'
######################
# HI_mid_avg
nc_HI_mid_avg = output_ncfile.createVariable('HI_mid_avg',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HI_mid_avg[:] = HI_mid_avg
# Add units and description attributes to the NetCDF variable if needed
nc_HI_mid_avg.units = 'haines index'
nc_HI_mid_avg.description = 'avg HI values - mid pressure level'
######################
# HI_high_avg
nc_HI_high_avg = output_ncfile.createVariable('HI_high_avg',np.float32, ('time','latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_HI_high_avg[:] = HI_high_avg
# Add units and description attributes to the NetCDF variable if needed
nc_HI_high_avg.units = 'haines index'
nc_HI_high_avg.description = 'avg HI values - high pressure level'
######################
# spatial_correlations_HIlow
nc_spatial_correlations_HIlow = output_ncfile.createVariable('spatial_correlations_HIlow',np.float32, ('latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_spatial_correlations_HIlow[:] = spatial_correlations_HIlow
# Add units and description attributes to the NetCDF variable if needed
nc_spatial_correlations_HIlow.units = '-1 1'
nc_spatial_correlations_HIlow.description = 'spatial correlation - prec - HI index'
######################
# spatial_correlations_HImid
nc_spatial_correlations_HImid = output_ncfile.createVariable('spatial_correlations_HImid',np.float32, ('latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_spatial_correlations_HImid[:] = spatial_correlations_HImid
# Add units and description attributes to the NetCDF variable if needed
nc_spatial_correlations_HImid.units = '-1 1'
nc_spatial_correlations_HImid.description = 'spatial correlation - prec - HI index'
######################
# spatial_correlations_HIhigh
nc_spatial_correlations_HIhigh = output_ncfile.createVariable('spatial_correlations_HIhigh',np.float32, ('latitude','longitude'))
# Assign data from 'tempC' to the NetCDF variable
nc_spatial_correlations_HIhigh[:] = spatial_correlations_HIhigh
# Add units and description attributes to the NetCDF variable if needed
nc_spatial_correlations_HIhigh.units = '-1 1'
nc_spatial_correlations_HIhigh.description = 'spatial correlation - prec - HI index'

# Close the NetCDF file
output_ncfile.close()