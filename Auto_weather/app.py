from flask import Flask, render_template, request, jsonify
import os
import requests
import traceback
import json  # Import the JSON module
from gfs_plot import generate_plot  # Import the plotting function
from datetime import datetime, timedelta  # Import datetime utilities

app = Flask(__name__)

# Load data from GFS.json
BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the current script
GFS_JSON_PATH = os.path.join(BASE_DIR, "GFS", "GFS.json")

try:
    with open(GFS_JSON_PATH, 'r') as f:
        data = json.load(f)
    print(f"Loaded variable data from {GFS_JSON_PATH}")
except Exception as e:
    print(f"Failed to load variable data from {GFS_JSON_PATH}: {e}")
    data = {}  # Fallback to an empty dictionary

BASE_URL = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl"

def get_latest_run():
    """
    Determine the most recent GFS run directory available.
    GFS runs are at 00, 06, 12, and 18 UTC. This function checks the current time
    and falls back to the previous run if the current one is not available.
    """
    now = datetime.utcnow()
    run_hours = [18, 12, 6, 0]  # GFS run hours in descending order
    for offset in range(2):  # Check today and yesterday
        check_date = now - timedelta(days=offset)
        date_str = check_date.strftime('%Y%m%d')  # Format as YYYYMMDD
        for run_hour in run_hours:
            run_dir = f"/gfs.{date_str}/{run_hour:02d}/atmos"
            # Construct a test URL to check if the run exists
            file_name = f"gfs.t{run_hour:02d}z.pgrb2.1p00.f000"
            test_url = f"{BASE_URL}?dir={run_dir}&file={file_name}&all_var=on&all_lev=on"
            try:
                response = requests.head(test_url, timeout=5)  # Use HEAD request to check availability
                if response.status_code == 200:
                    print(f"Found available run: {run_dir}")
                    return run_dir
                elif response.status_code == 403:
                    print(f"Run {run_dir} returned 403 Forbidden. Trying previous run.")
            except requests.exceptions.RequestException as e:
                print(f"Error checking run {run_dir}: {e}")
                continue  # Ignore errors and try the next run
    raise RuntimeError("No available GFS runs found in the last two days.")

@app.route('/')
def index():
    return render_template('index.html', variables=data)

@app.route('/plot', methods=['POST'])
def plot():
    # Get user input from the frontend
    variable = request.form.get('variable')
    time_step = request.form.get('time_step')
    color_scheme = request.form.get('color_scheme')  # Get the color scheme
    plot_type = request.form.get('plot_type')  # Get the plot type
    print(f"Received request to plot variable: {variable}, time step: {time_step}, color scheme: {color_scheme}, plot type: {plot_type}")

    # Handle custom colormap
    custom_colors = None
    if color_scheme == 'custom':
        color1 = request.form.get('color1')
        color2 = request.form.get('color2')
        color3 = request.form.get('color3')
        color4 = request.form.get('color4')
        custom_colors = [color1, color2, color3, color4]
        print(f"Custom colors selected: {custom_colors}")

    # Determine the typeOfLevel based on the variable
    type_of_level = None
    for level, variables in data.items():
        if variable in variables:
            type_of_level = level
            break

    if not type_of_level:
        error_message = f"Variable '{variable}' not found in the predefined data dictionary."
        print(error_message)
        return jsonify({'error': error_message})

    print(f"Determined typeOfLevel: {type_of_level} for variable: {variable}")

    # Determine the most recent available run directory
    try:
        run_dir = get_latest_run()
    except RuntimeError as e:
        error_message = f"Failed to determine the latest GFS run: {e}"
        print(error_message)
        return jsonify({'error': error_message})

    print(f"Using run directory: {run_dir}")

    # Extract the run hour from the run directory
    run_hour = run_dir.split('/')[-2]  # Extract the HH part from /gfs.YYYYMMDD/HH/atmos
    run_hour = run_hour[-2:]  # Get the last two characters (e.g., "00", "06", "12", "18")

    # Construct the GRIB2 file URL
    file_name = f"gfs.t{run_hour}z.pgrb2.1p00.f{time_step.zfill(3)}"
    params = {
        "dir": run_dir,
        "file": file_name,
        "all_var": "on",
        "all_lev": "on"
    }
    file_url = f"{BASE_URL}?dir={params['dir']}&file={params['file']}&all_var={params['all_var']}&all_lev={params['all_lev']}"
    print(f"Constructed GRIB2 file URL: {file_url}")

    # Ensure the static folder exists
    static_dir = os.path.join(BASE_DIR, 'static')
    if not os.path.exists(static_dir):
        os.makedirs(static_dir)
        print(f"Created static directory: {static_dir}")

    # Download the GRIB2 file
    grib_file_path = os.path.join(BASE_DIR, file_name)
    try:
        print(f"Downloading GRIB2 file to: {grib_file_path}")
        response = requests.get(file_url, stream=True)
        response.raise_for_status()
        with open(grib_file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Successfully downloaded GRIB2 file: {grib_file_path}")
    except requests.exceptions.RequestException as e:
        error_message = f"Failed to download GRIB2 file: {e}"
        print(error_message)
        return jsonify({'error': error_message})

    # Delegate the plotting to gfs_plot.py
    try:
        plot_path = generate_plot(
            grib_file_path, variable, type_of_level, time_step, static_dir, color_scheme, plot_type,
            custom_colors=custom_colors, run_dir=run_dir
        )
        print(f"Plot successfully generated at: {plot_path}")
    except Exception as e:
        error_message = f"Failed to generate plot: {e}\n{traceback.format_exc()}"
        print(error_message)
        return jsonify({'error': error_message})
    finally:
        # Delete the GRIB2 file after generating the plot
        if os.path.exists(grib_file_path):
            os.remove(grib_file_path)
            print(f"Deleted GRIB2 file: {grib_file_path}")

    return jsonify({'plot_url': f'/static/plot.png'})

if __name__ == '__main__':
    app.run(debug=True)
