#%%
from re import I
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import xarray as xr

xr_df = xr.open_dataset('data/gistemp1200_GHCNv4_ERSSTv5.nc')
xr_df

#%%
# Downsample the time series to yearly frequency
climate = xr_df.resample(time='Y').mean()
anomaly = climate['tempanomaly']

cbar_kwargs = {
    'orientation': 'horizontal',
    'fraction': 0.045,
    'pad': 0.01,
    'extend': 'neither'
}

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
ax.add_feature(NaturalEarthFeature('cultural', 'admin_0_countries', '10m'),
    facecolor='none', edgecolor='black')
ax.set_extent([-150, 150, -55, 85])

i = -1
date = pd.to_datetime(anomaly.isel(time=i)['time'].values)
ax.set_title("Temperature anomaly in " + str(date.year) + " [°C]")
anomaly.isel(time=i).plot.imshow(ax=ax, add_labels=False, add_colorbar=True,
    vmin=-4, vmax=4, cmap='coolwarm',
    cbar_kwargs=cbar_kwargs, interpolation='bicubic')
plt.savefig("global_map.png", bbox_inches='tight', dpi=150)
plt.show

