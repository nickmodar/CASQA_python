################################################################################
################################################################################
################## ~ Season 1: Haines Index Figure (low) ~ #####################
################################################################################
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import geopandas as gpd
import cartopy.crs as ccrs
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import griddata
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

fn1 = 'C:/Users/modarn/Documents/CASQA/data/ta_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc'
season2_Fin = 'C:/Users/modarn/Documents/CASQA/data/season1997_1998.nc'
season2_in = nc.Dataset(season2_Fin)
data1IN = nc.Dataset(fn1)# load in original nc file for lat,lon (and time maybe) dimensions.
lat = data1IN['lat'][:]
lon = data1IN['lon'][:]
HI_low_in = season2_in['HI_low'][:]
HI_low = np.mean(HI_low_in, axis=0)

# Load California boundaries from GeoJSON file
california = gpd.read_file('C:/Users/modarn/Documents/CASQA/CASQA_python/output.geojson')

# Define the target projection (EPSG 4326 - WGS 84)
target_projection = ccrs.PlateCarree()

# Get the extent of California boundaries
california_extent = california.total_bounds

# Create a figure and axis with Cartopy projection, setting the extent to California's boundaries
padding = 1.0
fig, ax = plt.subplots(figsize=(11, 9), subplot_kw={'projection': target_projection}, constrained_layout=True)
ax.set_extent([california_extent[0], california_extent[2], california_extent[1], california_extent[3]], crs=target_projection)

# Plot the California boundaries with black border
ax.add_geometries(california['geometry'], crs=target_projection, edgecolor='black', facecolor='none', linewidth=1)

# Define the meshgrid for latitudes and longitudes within California's extent
lon_mesh, lat_mesh = np.meshgrid(np.linspace(california_extent[0], california_extent[2], 192),
                                 np.linspace(california_extent[1], california_extent[3], 96))

# Create a mask to hide data outside of California's outline
mask = california.geometry.unary_union

# Apply the mask to corrHI_low , np.nan
HI_low_masked = np.where(mask.contains(gpd.points_from_xy(lon_mesh.flatten(), lat_mesh.flatten())), HI_low.flatten(), np.nan)
HI_low_masked = HI_low_masked.reshape(lon_mesh.shape)

# Create a grid of all coordinates within California
lon_grid, lat_grid = np.meshgrid(
    np.linspace(california_extent[0], california_extent[2], 192),
    np.linspace(california_extent[1], california_extent[3], 96))

# Create contour plot for the masked data
cmap = plt.get_cmap('RdYlGn')
contour = plt.contourf(lon_mesh, lat_mesh, HI_low_masked,
                       levels= 10, extend='both', cmap=cmap.reversed('RdYlGn'))#, vmin=0, vmax=6)
#cbar = plt.colorbar(contour, pad=0.03, orientation='vertical')
# Set custom tick locations and labels
#tick_labels = ['0','1','2','3','4','5']
#cbar.set_ticklabels(tick_labels)
#cbar.set_label('Haines Index')

# Set the extent to match California boundaries
ax.set_extent([california_extent[0] - padding, california_extent[2] + padding,
               california_extent[1] - padding, california_extent[3] + padding], crs=target_projection)

# Add gridlines with labels
#gl = ax.gridlines(crs=target_projection, draw_labels=True, linewidth=0.5, linestyle='--', color='gray')
#gl.xlabels_top = False
#gl.ylabels_right = False
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER

# Add a title
#ax.set_title('Correlation Data within California Boundaries (Contour Plot)')

# Show the plot
plt.show()

