import cfgrib
import json

def read_variables(file_path, output_json):
    """
    Reads a GRIB2 file and saves all available variables along with their units into a JSON file, resolving multiple values for unique keys.

    :param file_path: The path to the GRIB2 file.
    :param output_json: The path to the output JSON file.
    """
    type_of_level_filters = [
        'meanSea', 'hybrid', 'atmosphere', 'surface', 'planetaryBoundaryLayer',
        'isobaricInPa', 'isobaricInhPa', 'heightAboveGround', 'depthBelowLandLayer',
        'heightAboveSea', 'atmosphereSingleLayer', 'lowCloudLayer', 'middleCloudLayer',
        'highCloudLayer', 'cloudCeiling', 'heightAboveGroundLayer', 'tropopause',
        'maxWind', 'isothermZero', 'highestTroposphericFreezing', 'pressureFromGroundLayer',
        'sigmaLayer', 'sigma', 'potentialVorticity'
    ]

    variables_by_level = {}

    for level in type_of_level_filters:
        try:
            ds = cfgrib.open_dataset(file_path, filter_by_keys={'typeOfLevel': level})
            # Filter variables to include only those after 'valid_time'
            variables = list(ds.variables.keys())
            if 'valid_time' in variables:
                valid_time_index = variables.index('valid_time')
                variables = variables[valid_time_index + 1:]
            else:
                variables = []

            # Collect variables with their units
            variables_with_units = {}
            for var in variables:
                var_obj = ds.variables[var]
                units = var_obj.attrs.get('units', 'unknown')  # Get units or default to 'unknown'
                variables_with_units[var] = units

            variables_by_level[level] = variables_with_units
        except Exception as e:
            variables_by_level[level] = f"Error: {str(e)}"

    # Save the variables with units to a JSON file
    with open(output_json, 'w') as json_file:
        json.dump(variables_by_level, json_file, indent=4)

    print(f"Variables with units saved to {output_json}")

if __name__ == "__main__":
    # Path to the GRIB2 file
    file_path = "gfs.t18z.pgrb2.1p00.f000"

    # Path to the output JSON file
    output_json = "variables_with_units.json"

    # Read and save variables with units
    read_variables(file_path, output_json)
