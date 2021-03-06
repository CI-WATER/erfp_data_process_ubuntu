# erfp_data_process_ubuntu
Code to use to prepare input data for RAPID from ECMWF forecast on Linux (created in Ubuntu) using Multiprocessing.
##Step 1: Install RAPID
**For Ubuntu:**
```
$ apt-get install gfortran g++
```
Follow the instructions on page 10-14: http://rapid-hub.org/docs/RAPID_Azure.pdf.

Add run_rapid.sh to the rapid/run directory with the lines:
```
cd /home/alan/work/rapid/run/
./rapid
```

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
##Step 3: Install tethys_dataset_services
```
$ pip install requests_toolbelt
$ pip install tethys_dataset_services
```
##Step 4: Download the source code
```
$ cd /path/to/your/scripts/
$ git clone https://github.com/CI-WATER/erfp_data_process_ubuntu.git
$ git submodule init
$ git submodule update
```
##Step 5: Create folders for RAPID input and for downloading ECMWF
In this instance:
```
$ cd /home/alan/
$ mkdir work/rapid/input work/ecmwf work/logs
```
##Step 6: Change the locations in the files
Go into *rapid_process_async_ubuntu.py* and change these variables for your instance:
```python
    rapid_files_location = '/home/alan/work/rapid'
    ecmwf_forecast_location = "/home/alan/work/ecmwf"
    ckan_api_endpoint = 'http://ciwckan.chpc.utah.edu'
    ckan_api_key = 'areally-good-key'
```
Go into *rapid_process.sh* and change make sure the path locations and variables are correct for your instance.
##Step 7: Make sure permissions are correct for these files and any directories the script will use

Example:
```
$ chmod 554 rapid_process_async_ubuntu.py
$ chmod 554 rapid_process.sh
```
##Step 8: Add RAPID files to the work/rapid/input directory
Example:
```
$ ls work/rapid/input
huc_region_1209
$ ls -lh work/rapid/input/huc_region_1209
-r--r--r-- 1 alan alan 163K Mar  6 10:01 k.csv
-r--r--r-- 1 alan alan 163K Mar  6 10:01 kfac.csv
-r--r--r-- 1 alan alan 340K Mar  6 19:22 rapid_connect.csv
-rw-r--r-- 1 alan alan 5.1K Mar 25 04:15 rapid_namelist_huc_4_1209.dat
-r--r--r-- 1 alan alan  99K Mar  9 07:52 riv_bas_id_huc_4_1209.csv
-r--r--r-- 1 alan alan 1.5M Mar  9 08:03 weight_high_res.csv
-r--r--r-- 1 alan alan 1.2M Mar  9 08:03 weight_low_res.csv
-r--r--r-- 1 alan alan  55K Mar  6 10:01 x.csv
```
##Step 9: Create CRON job to run the scripts twice daily
See: http://askubuntu.com/questions/2368/how-do-i-set-up-a-cron-job

You only need to run rapid_process.sh
```
$ ./rapid_process.sh
```
#Troubleshooting
If you see this error:
ImportError: No module named packages.urllib3.poolmanager
```
$ pip install pip --upgrade
```
Restart your terminal
```
$ pip install requests --upgrade
```
