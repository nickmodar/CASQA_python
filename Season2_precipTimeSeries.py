import netCDF4 as nc
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
#Convert from kg/m^2/s to kg/m^2/day.
#Convert from kg/m^2/day to inches/day.
#Convert from kg/m^2/s to kg/m^2/day:
#Since there are 86,400 seconds in a day (24 hours * 60 minutes * 60 seconds), you can multiply the precipitation flux in kg/m^2/s by 86,400 to get kg/m^2/day.
#Convert from kg/m^2/day to inches/day:
#To convert from kg/m^2 to inches, you can use the following conversion factor:
#1 kg/m^2 ≈ 0.0393701 inches
#So, multiply the precipitation flux in kg/m^2/day by 0.0393701 to get inches/day.
#Here's the formula for the conversion:
#Precipitation_flux (inches/day) = Precipitation_flux (kg/m^2/s) * 86,400 (s/day) * 0.0393701 (inches/kg/m^2)
####################################################################################################
# LOAD IN NETCDF FILES
fn3 = 'C:/Users/modarn/Documents/CASQA/data/pr_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc'
data3IN = nc.Dataset(fn3)
## define varaibles from netCDF files
prec_in = data3IN['pr'][:]
conversion_factor = 3401.2
prec_inches_per_day = prec_in * conversion_factor
lat = data3IN['lat'][:]
lon = data3IN['lon'][:]
## TIME - SEASON 1
# 1982/1983 season, fire season = May through August
# Extract the time data and convert it to datetime objects
time_data = data3IN.variables['time'][:]
time_units = data3IN.variables['time'].units
time_calendar = data3IN.variables['time'].calendar
time_datetime = nc.num2date(time_data, units=time_units, calendar=time_calendar)
# Wet season - define time
# Define the date range from November 1, 1997 to March 31, 1998
start_date_wet = datetime(1997, 11, 1)
end_date_wet = datetime(1998, 3, 31)
# Find the indices of time data that fall within the specified date range
time_indices_wet = [idx for idx, t in enumerate(time_datetime) if start_date_wet <= t <= end_date_wet]
# Create a variable that shows the datetime values selected by time_indices to check work
selected_datetime_values_wet = [time_datetime[idx] for idx in time_indices_wet]
# Extract the data for the specified date range for PRECIP, TEMPERATURE AND RELATIVE HUMIDITY VARIABLE
prec = prec_inches_per_day[time_indices_wet, :, :]
#####################################################################
# Define the latitude and longitude boundaries for California
lat_min, lat_max = 32.5342307609976, 42.00965914828148
lon_min, lon_max = 114.13445790587905, 124.41060660766607
# Calculate the corresponding indices for the lat and lon slices
lat_indices = np.where((data3IN['lat'][:] >= lat_min) & (data3IN['lat'][:] <= lat_max))[0]
lon_indices = np.where((data3IN['lon'][:] >= lon_min) & (data3IN['lon'][:] <= lon_max))[0]

# Subset the precip_data based on California's boundaries
california_precip_data = prec[:, lat_indices, lon_indices]
# Create an empty list to store the summed precipitation values
daily_sum_precip = []
# Loop through each row in california_precip_data
for day_data in california_precip_data:
    # Sum the 6 data points for each day and append to the list
    daily_sum = np.ma.sum(day_data)
    daily_sum_precip.append(daily_sum)
# Convert the list to a NumPy array
daily_sum_precip = np.array(daily_sum_precip)
# Create an empty list to store the summed precipitation values
daily_sum_precip = []
# Loop through each row in california_precip_data
for day_data in california_precip_data:
    # Sum the 6 data points for each day and append to the list
    daily_sum = np.ma.sum(day_data)
    daily_sum_precip.append(daily_sum)
# Convert the list to a NumPy array
daily_sum_precip = np.array(daily_sum_precip)
# Add two additional data points to daily_sum_precip
additional_data_points = [0.0, 0.0]  # Add your desired values here
daily_sum_precip = np.concatenate((daily_sum_precip, additional_data_points))

#################################

# DATETIME FOR PLOT
start_date_wet = datetime(1997, 11, 1)
end_date_wet = datetime(1998, 3, 31)
dates_list = [start_date_wet + timedelta(days=x) for x in range(0, (end_date_wet-start_date_wet).days+2)]
myFmt = mdates.DateFormatter('%b %Y')
months = mdates.MonthLocator()
days = mdates.DayLocator(bymonthday=(1))
dayss = mdates.DayLocator()

# PLOT THE FIGURE
plt.figure(figsize=(12, 8))
plt.subplots_adjust(bottom=0.2)
ax = plt.gca()
ax.plot(dates_list, daily_sum_precip, linestyle='-',linewidth='5', color='royalblue')

# Y-AXIS
min_value = np.min(daily_sum_precip)
max_value = np.max(daily_sum_precip)
y_ticks = np.linspace(0, 5,  6)
ax.set_ylim(0, 5)
plt.yticks(y_ticks)

# X-AXIS
ax.set_xlim(dates_list[0], dates_list[-1])
ax.set_xlim(dates_list[0],dates_list[-1])
ax.xaxis.set_major_formatter(myFmt)
ax.xaxis.set_major_locator(days)
ax.xaxis.set_minor_locator(dayss)
locs = list(ax.get_xticks())+ [mdates.date2num(datetime(1998,4,1))]
locator= matplotlib.ticker.FixedLocator(locs)
ax.xaxis.set_major_locator(locator)

# PLOT ATTRIBUTES
# Increase font size for tick labels, axis labels, and title
ax.tick_params(axis='both', which='both', labelsize=16)  # Adjust the labelsize as needed
ax.set_xlabel(r'Event 2: California Wet Season' + '\n' + r'November 1, 1997 - April 1, 1998', fontsize=24)
ax.set_ylabel('Average Precipitation'+ '\n' +  '(inches/day)', fontsize=24)  # Adjust the fontsize as needed

plt.show()