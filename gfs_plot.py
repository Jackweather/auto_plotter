import matplotlib.pyplot as plt
import cfgrib
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.colors as mcolors  # Import Normalize for colormap normalization
from scipy.interpolate import griddata


# Define unit conversion mapping for USA units with unit labels
UNIT_CONVERSIONS = {
    "meanSea": {
        "prmsl": {"func": lambda x: x / 100, "unit": "hPa"},  # Pa to hPa
        "mslet": {"func": lambda x: x / 100, "unit": "hPa"},  # Pa to hPa
    },
    "hybrid": {
        "clwmr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "icmr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "rwmr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "snmr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "grle": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
    },
    "atmosphere": {
        "refc": {"func": lambda x: x, "unit": "dB"},  # dB (no conversion needed)
        "tcc": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
    },
    "surface": {
        "vis": {"func": lambda x: x / 1609.34, "unit": "miles"},  # m to miles
        "gust": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "hindex": {"func": lambda x: x, "unit": ""},  # Numeric (no conversion needed)
        "sp": {"func": lambda x: x / 100, "unit": "hPa"},  # Pa to hPa
        "orog": {"func": lambda x: x * 3.28084, "unit": "feet"},  # m to feet
        "t": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
        "cnwat": {"func": lambda x: x, "unit": "kg/m²"},  # kg/m² (no conversion needed)
        "sdwe": {"func": lambda x: x, "unit": "kg/m²"},  # kg/m² (no conversion needed)
        "sde": {"func": lambda x: x * 3.28084, "unit": "feet"},  # m to feet
        "sithick": {"func": lambda x: x * 3.28084, "unit": "feet"},  # m to feet
        "cpofp": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
        "prate": {"func": lambda x: x * 3600, "unit": "mm/hr"},  # kg/m²/s to mm/hr
        "fsr": {"func": lambda x: x, "unit": "m"},  # m (no conversion needed)
        "fricv": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "veg": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
        "wilt": {"func": lambda x: x, "unit": ""},  # Fraction (no conversion needed)
        "fldcp": {"func": lambda x: x, "unit": ""},  # Fraction (no conversion needed)
        "SUNSD": {"func": lambda x: x, "unit": "s"},  # s (no conversion needed)
        "lftx": {"func": lambda x: x, "unit": "K"},  # K (no conversion needed)
        "cape": {"func": lambda x: x, "unit": "J/kg"},  # J/kg (no conversion needed)
        "cin": {"func": lambda x: x, "unit": "J/kg"},  # J/kg (no conversion needed)
        "lftx4": {"func": lambda x: x, "unit": "K"},  # K (no conversion needed)
        "lsm": {"func": lambda x: x, "unit": ""},  # (0-1) (no conversion needed)
        "siconc": {"func": lambda x: x, "unit": ""},  # (0-1) (no conversion needed)
        "sit": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
    },
    "planetaryBoundaryLayer": {
        "u": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "v": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "VRATE": {"func": lambda x: x, "unit": "m²/s"},  # m²/s (no conversion needed)
    },
    "isobaricInPa": {
        "gh": {"func": lambda x: x * 0.0328084, "unit": "feet"},  # gpm to feet
        "t": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
        "r": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
        "q": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "w": {"func": lambda x: x, "unit": "Pa/s"},  # Pa/s (no conversion needed)
        "wz": {"func": lambda x: x, "unit": "m/s"},  # m/s (no conversion needed)
        "u": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "v": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "absv": {"func": lambda x: x, "unit": "s⁻¹"},  # s⁻¹ (no conversion needed)
        "o3mr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
    },
    "isobaricInhPa": {
        "gh": {"func": lambda x: x * 0.0328084, "unit": "feet"},  # gpm to feet
        "t": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
        "r": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
        "q": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
        "w": {"func": lambda x: x, "unit": "Pa/s"},  # Pa/s (no conversion needed)
        "wz": {"func": lambda x: x, "unit": "m/s"},  # m/s (no conversion needed)
        "u": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "v": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "absv": {"func": lambda x: x, "unit": "s⁻¹"},  # s⁻¹ (no conversion needed)
        "o3mr": {"func": lambda x: x, "unit": "kg/kg"},  # kg/kg (no conversion needed)
    },
    "heightAboveGround": {
        "refd": {"func": lambda x: x, "unit": "dB"},  # dB (no conversion needed)
    },
    "depthBelowLandLayer": {
        "st": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
        "soilw": {"func": lambda x: x, "unit": ""},  # Proportion (no conversion needed)
        "soill": {"func": lambda x: x, "unit": ""},  # Proportion (no conversion needed)
    },
    "atmosphereSingleLayer": {
        "pwat": {"func": lambda x: x, "unit": "kg/m²"},  # kg/m² (no conversion needed)
        "cwat": {"func": lambda x: x, "unit": "kg/m²"},  # kg/m² (no conversion needed)
        "r": {"func": lambda x: x, "unit": "%"},  # % (no conversion needed)
        "tozne": {"func": lambda x: x, "unit": "DU"},  # DU (no conversion needed)
    },
    "tropopause": {
        "trpp": {"func": lambda x: x / 100, "unit": "hPa"},  # Pa to hPa
        "icaht": {"func": lambda x: x * 3.28084, "unit": "feet"},  # m to feet
        "gh": {"func": lambda x: x * 0.0328084, "unit": "feet"},  # gpm to feet
        "t": {"func": lambda x: x - 273.15, "unit": "°C"},  # K to °C
        "u": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "v": {"func": lambda x: x * 2.23694, "unit": "mph"},  # m/s to mph
        "vwsh": {"func": lambda x: x, "unit": "s⁻¹"},  # s⁻¹ (no conversion needed)
    },
    # Add other levels and variables as needed...
}

def convert_units(data, variable, type_of_level):
    """
    Convert data to USA units based on the variable and type_of_level.

    :param data: The data array to convert.
    :param variable: The variable name.
    :param type_of_level: The type of level.
    :return: Converted data.
    """
    if type_of_level in UNIT_CONVERSIONS and variable in UNIT_CONVERSIONS[type_of_level]:
        conversion_func = UNIT_CONVERSIONS[type_of_level][variable]["func"]
        return conversion_func(data)  # The conversion function is applied to the data
    return data  # Return unmodified data if no conversion is defined

def get_unit_label(variable, type_of_level):
    """
    Retrieve the unit label for the given variable and type_of_level.

    :param variable: The variable name.
    :param type_of_level: The type of level.
    :return: The unit label as a string.
    """
    if type_of_level in UNIT_CONVERSIONS and variable in UNIT_CONVERSIONS[type_of_level]:
        return UNIT_CONVERSIONS[type_of_level][variable]["unit"]
    return "Unknown Unit"

def generate_plot(grib_file_path, variable, type_of_level, time_step, static_dir, color_scheme, plot_type, step_type="instant", custom_colors=None, run_dir=None):

    print(f"Opening GRIB2 file: {grib_file_path}")
    try:
        ds = cfgrib.open_dataset(
            grib_file_path,
            filter_by_keys={'typeOfLevel': type_of_level, 'stepType': step_type}
        )
    except cfgrib.dataset.DatasetBuildError as e:
        if "multiple values for unique key" in str(e):
            raise ValueError(
                "The GRIB2 file contains multiple values for 'stepType'. "
                "Please specify one of the following step types: 'instant', 'avg', or 'accum'."
            ) from e
        else:
            raise

    if variable not in ds:
        raise ValueError(f"Variable '{variable}' not found. Available: {list(ds.data_vars)}")

    # Load the data
    data = ds[variable].values
    lats = ds.latitude.values
    lons = ds.longitude.values

    # Debugging: Print raw data range
    print(f"Raw Data Min: {data.min()}, Max: {data.max()}")

    # Convert data to USA units
    data = convert_units(data, variable, type_of_level)

    # Debugging: Print converted data range
    print(f"Converted Data Min: {data.min()}, Max: {data.max()}")

    print("Loading county boundaries...")
    counties_path = os.path.join(os.path.dirname(__file__), "counties.json")
    counties = gpd.read_file(counties_path)

    print("Rendering map...")
    fig = plt.figure(figsize=(14, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.4)
    ax.add_feature(cfeature.LAKES, edgecolor="black", linewidth=0.3)
    ax.add_feature(cfeature.RIVERS, linewidth=0.3)

    # Define the map extent
    extent = [-130, -60, 20, 55]  # [west, east, south, north]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Draw counties (FAST)
    ax.add_geometries(
        counties.geometry,
        crs=ccrs.PlateCarree(),
        facecolor="none",
        edgecolor="gray",
        linewidth=0.3
    )

    # Debugging: Print data shape and range
    print(f"Data shape: {data.shape}, Min: {data.min()}, Max: {data.max()}")
    print(f"Lats shape: {lats.shape}, Lons shape: {lons.shape}")

    # Fix longitudes if they are in the range 0–360
    if lons.max() > 180:
        lons = np.where(lons > 180, lons - 360, lons)

    # Apply Gaussian smoothing to reduce noise
    data = gaussian_filter(data, sigma=1)

    # Mask invalid or missing data
    data = np.ma.masked_invalid(data)

    # Check if the data is entirely NaN
    if np.all(np.isnan(data)):
        raise ValueError(
            f"Data for variable '{variable}' at level '{type_of_level}' contains only NaN values. "
            "Cannot generate plot."
        )

    # Ensure longitude and latitude grid alignment
    if lats.ndim == 1 and lons.ndim == 1:
        lons, lats = np.meshgrid(lons, lats)  # Create 2D mesh grids

    # Interpolate data to a higher resolution grid
    grid_x, grid_y = np.linspace(lons.min(), lons.max(), 500), np.linspace(lats.min(), lats.max(), 500)
    grid_lons, grid_lats = np.meshgrid(grid_x, grid_y)

    # Flatten the input arrays for griddata
    points = np.array([lons.flatten(), lats.flatten()]).T
    values = data.flatten()
    data = griddata(points, values, (grid_lons, grid_lats), method='linear')

    # Update lons and lats to the new grid
    lons, lats = grid_lons, grid_lats

    # Handle NaN values in the interpolated data
    if np.any(np.isnan(data)):
        print("Warning: Interpolated data contains NaN values. Replacing NaN values with the minimum valid value.")
        data = np.nan_to_num(data, nan=np.nanmin(data))  # Replace NaN with the minimum valid value

    # Define contour levels
    data_min, data_max = data.min(), data.max()
    levels = np.linspace(data_min, data_max, 50)

    # Define custom colormaps
    custom_colormaps = {
        "white_to_blue": mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["white", "blue"]),
        "white_to_red": mcolors.LinearSegmentedColormap.from_list("white_to_red", ["white", "red"]),
        "white_to_green": mcolors.LinearSegmentedColormap.from_list("white_to_green", ["white", "green"]),
        "white_to_purple": mcolors.LinearSegmentedColormap.from_list("white_to_purple", ["white", "purple"]),
        "white_to_orange": mcolors.LinearSegmentedColormap.from_list("white_to_orange", ["white", "orange"]),
    }

    # Handle custom colormap
    if color_scheme == 'custom' and custom_colors:
        cmap = mcolors.LinearSegmentedColormap.from_list("custom", custom_colors)
    elif color_scheme in custom_colormaps:
        cmap = custom_colormaps[color_scheme]
    else:
        raise ValueError(f"Invalid color scheme: {color_scheme}")

    # Define normalization for the colormap
    norm = mcolors.Normalize(vmin=data_min, vmax=data_max)

    # Use the selected plot type
    if plot_type == "contourf":
        plot = ax.contourf(
            lons,
            lats,
            data,
            levels=levels,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree()
        )
    elif plot_type == "contour":
        contour = ax.contour(
            lons,
            lats,
            data,
            levels=levels,
            colors="black",
            linewidths=0.5,
            transform=ccrs.PlateCarree()
        )
        # Add labels to all contour lines
        ax.clabel(
            contour,
            inline=True,
            inline_spacing=5,  # Adjust spacing for better label placement
            fontsize=8,
            fmt="%.1f"  # Format labels to one decimal place
        )
    elif plot_type == "pcolormesh":
        plot = ax.pcolormesh(
            lons,
            lats,
            data,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree()
        )
    else:
        raise ValueError(f"Invalid plot type: {plot_type}")

    # Add a colorbar for contourf and pcolormesh
    if plot_type in ["contourf", "pcolormesh"]:
        cbar = plt.colorbar(plot, orientation="horizontal", pad=0.03, aspect=50)
        # Add the unit label to the color bar
        unit_label = get_unit_label(variable, type_of_level)
        cbar.set_label(f"{variable} ({unit_label})", fontsize=10)

    # Extract yymmdd and run hour from run_dir
    if run_dir:
        run_date = run_dir.split('/')[1][4:]  # Extract yymmdd from /gfs.YYYYMMDD/HH/atmos
        run_hour = run_dir.split('/')[2]  # Extract hh from /gfs.YYYYMMDD/HH/atmos
        run_info = f"Run: {run_date} {run_hour}Z"
    else:
        run_info = "Run: Unknown"

    # Set the title for the plot
    ax.set_title(f"{variable}   •   Forecast Hour {time_step}\n{run_info}", fontsize=14, fontweight="bold")

    # Save the plot
    os.makedirs(static_dir, exist_ok=True)
    plot_path = os.path.join(static_dir, "plot.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")  # Increase DPI for better quality
    plt.close()

    print("Saved:", plot_path)
    return plot_path
