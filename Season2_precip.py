################################################################################
################################################################################
################## ~ Season 1: Average Precipitation  ~ ########################
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
prec_in = season2_in['precipitation'][:]
prec = np.mean(prec_in, axis=0)

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
prec_masked = np.where(mask.contains(gpd.points_from_xy(lon_mesh.flatten(), lat_mesh.flatten())), prec.flatten(), np.nan)
prec_masked = prec_masked.reshape(lon_mesh.shape)

# Create a grid of all coordinates within California
lon_grid, lat_grid = np.meshgrid(
    np.linspace(california_extent[0], california_extent[2], 192),
    np.linspace(california_extent[1], california_extent[3], 96))

# Create contour plot for the masked data

nws_precip_colors = [
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6",# 8.00 - 10.00 inches
    "#fdfdfd"   # 10.00+
]
precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

# Create a contour plot with the custom colormap
contour = plt.contourf(lon_mesh, lat_mesh, prec_masked,
                       levels=10, extend='both', cmap=precip_colormap)
cbar = plt.colorbar(contour, pad=0.03, orientation='vertical')
#color_levels = [0.01, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
#color_labels = [str(level) for level in color_levels]
# Set custom tick locations and labels
#tick_locations = np.arange(len(color_levels))
#cbar.set_ticks(tick_locations)
#cbar.set_ticklabels(color_labels)
# Set color bar label
#cbar.set_label('Precipitation (inches)')


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

