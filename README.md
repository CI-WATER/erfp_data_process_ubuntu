# erfp_data_process_ubuntu
Code to use to prepare input data for RAPID from ECMWF forecast on Linux (created in Ubuntu) using Multiprocessing.
##Step 1: Install RAPID
**For Ubuntu:**
```
$ apt-get install gfortran g++
```
Follow the instructions on page 10-14: http://rapid-hub.org/docs/RAPID_Azure.pdf.

##Step 2: Install netCDF4-python
###Install on Ubuntu:
```
$ apt-get install python-dev zlib1g-dev libhdf5-serial-dev libnetcdf-dev 
$ pip install numpy
$ pip install netCDF4
```
###Install on Redhat:
*Note: this tool was desgined and tested in Ubuntu*
```
$ yum install netcdf4-python
$ yum install hdf5-devel
$ yum install netcdf-devel
$ pip install numpy
$ pip install netCDF4
```
##Step 3: Download the source code
```
$ cd /path/to/your/scripts/
$ git clone https://github.com/CI-WATER/erfp_data_process_ubuntu.git
```
##Step 4: Change the locations in the files
Go into *rapid_process_async_ubuntu.py* and change these variables to the appropriate locations:
```
    rapid_files_location = '/home/alan/work/rapid'
    ecmwf_forecast_location = "/home/alan/work/ecmwf"
    ckan_api_endpoint = 'http://ciwckan.chpc.utah.edu'
    ckan_api_key = 'areally-good-key'
```
Go into *rapid_process.sh* and change make sure the path locations and variables are correct for your instance.
##Step 5: Make sure permissions are correct for these files and any directories the script will use

Example:
```
chmod 554 rapid_process_async_ubuntu.py
chmod 554 rapid_process.sh
```
##Step 6: Add RAPID files to the work/rapid/input directory
Example:
```
$ ls work/rapid/input
huc_region_1209
$ ls work/rapid/input/huc_region_1209
-rwxr-xr-x 1 alan alan 163K Mar  6 10:01 k.csv
-rwxr-xr-x 1 alan alan 163K Mar  6 10:01 kfac.csv
-rwxr-xr-x 1 alan alan 340K Mar  6 19:22 rapid_connect.csv
-rw------- 1 alan alan 5.1K Mar 25 04:15 rapid_namelist_huc_4_1209.dat
-rwxr-xr-x 1 alan alan  99K Mar  9 07:52 riv_bas_id_huc_4_1209.csv
-rw-r--r-- 1 alan alan 1.5M Mar  9 08:03 weight_high_res.csv
-rw-r--r-- 1 alan alan 1.2M Mar  9 08:03 weight_low_res.csv
-rwxr-xr-x 1 alan alan  55K Mar  6 10:01 x.csv
```
##Step 7: Create CRON job to run the scripts twice daily
See: http://askubuntu.com/questions/2368/how-do-i-set-up-a-cron-job

You only need to run rapid_process.sh
```
./rapid_process.sh
```
