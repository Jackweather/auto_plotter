import requests

def download_file(url, output_path):
    """
    Downloads a file from the given URL and saves it to the specified output path.

    :param url: The URL of the file to download.
    :param output_path: The local path where the file will be saved.
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an error for HTTP issues

        with open(output_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        print(f"File downloaded successfully and saved to {output_path}")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while downloading the file: {e}")

if __name__ == "__main__":
    # URL of the file to download
    url = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl?dir=%2Fgfs.20260107%2F18%2Fatmos&file=gfs.t18z.pgrb2.1p00.f000&all_var=on&all_lev=on"
    
    # Local path to save the downloaded file
    output_path = "gfs.t18z.pgrb2.1p00.f000"

    # Download the file
    download_file(url, output_path)
