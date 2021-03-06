import csv
import datetime
from glob import glob
import itertools
import multiprocessing
import netCDF4 as NET
import numpy as np
import os
import subprocess
from tempfile import mkstemp
from time import sleep
from shutil import move, rmtree
import ftp_ecmwf_download
from CreateInflowFileFromECMWFRunoff import CreateInflowFileFromECMWFRunoff
from sfpt_dataset_manager.dataset_manager import ECMWFRAPIDDatasetManager

#------------------------------------------------------------------------------
#functions
#------------------------------------------------------------------------------
def find_current_rapid_output(forecast_directory, basin_name):
    """
    Finds the most current files output from RAPID
    """
    if os.path.exists(forecast_directory):
        basin_files = glob(os.path.join(forecast_directory,"Qout_%s_*.nc" % basin_name))
        if len(basin_files) >0:
            return basin_files
    #there are none found
    return None

def csv_to_list(csv_file, delimiter=','):
    """
    Reads in a CSV file and returns the contents as list,
    where every row is stored as a sublist, and each element
    in the sublist represents 1 cell in the table.

    """
    with open(csv_file, 'rb') as csv_con:
        reader = csv.reader(csv_con, delimiter=delimiter)
        return list(reader)

def get_comids_in_netcdf_file(reach_id_list, prediction_file):
    """
    Gets the subset comid_index_list, reordered_comid_list from the netcdf file
    """
    data_nc = NET.Dataset(prediction_file, mode="r")
    com_ids = data_nc.variables['COMID'][:]
    data_nc.close()
    try:
        #get where comids are in netcdf file
        netcdf_reach_indices_list = np.where(np.in1d(com_ids, reach_id_list))[0]
    except Exception as ex:
        print ex

    return netcdf_reach_indices_list, com_ids[netcdf_reach_indices_list]

def compute_initial_rapid_flows(prediction_files, basin_name, input_directory, forecast_date_timestep):
    """
    Gets mean of all 52 ensembles 12-hrs in future and prints to csv as initial flow
    Qinit_file (BS_opt_Qinit)
    The assumptions are that Qinit_file is ordered the same way as rapid_connect_file
    if subset of list, add zero where there is no flow
    """
    #remove old init files for this basin
    past_init_flow_files = glob(os.path.join(input_directory, 'Qinit_file_%s_*.csv' % basin_name))
    for past_init_flow_file in past_init_flow_files:
        try:
            os.remove(past_init_flow_file)
        except:
            pass

    init_file_location = os.path.join(input_directory,'Qinit_file_%s_%s.csv' % (basin_name, forecast_date_timestep))
    #check to see if exists and only perform operation once
    if prediction_files:
        #get list of COMIDS
        connectivity_file = csv_to_list(os.path.join(input_directory,'rapid_connect.csv'))
        comid_list = np.array([int(row[0]) for row in connectivity_file])


        print "Finding COMID indices ..."
        comid_index_list, reordered_comid_list = get_comids_in_netcdf_file(comid_list, prediction_files[0])
        print "Extracting Data ..."
        reach_prediciton_array = np.zeros((len(comid_list),len(prediction_files),1))
        #get information from datasets
        for file_index, prediction_file in enumerate(prediction_files):
            try:
                #Get hydrograph data from ECMWF Ensemble
                data_nc = NET.Dataset(prediction_file, mode="r")
                qout_dimensions = data_nc.variables['Qout'].dimensions
                if qout_dimensions[0] == 'Time' and qout_dimensions[1] == 'COMID':
                    data_values_2d_array = data_nc.variables['Qout'][2,comid_index_list]
                    for comid_index, comid in enumerate(reordered_comid_list):
                        reach_prediciton_array[comid_index][file_index] = data_values_2d_array[comid_index]
                else:
                    print "Invalid ECMWF forecast file", prediction_file
                    data_nc.close()
                    continue
                data_nc.close()

            except Exception, e:
                print e
                #pass

        print "Analyzing Data ..."
        output_data = []
        for comid in comid_list:
            try:
                #get where comids are in netcdf file
                comid_index = np.where(reordered_comid_list==comid)[0][0]
            except Exception:
                #comid not found in list. Adding zero init flow ...
                output_data.append([0])
                pass
                continue

            #get mean of series as init flow
            output_data.append([np.mean(reach_prediciton_array[comid_index])])

        print "Writing Output ..."
        with open(init_file_location, 'wb') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(output_data)
    else:
        print "No ensembles for basin found"
     
     
def update_namelist_file(namelist_file, output_directory, vlat_file, 
    qout_file, qinit_file, ensemble_number):
    """
    Update RAPID namelist file with new inflow file and output location
    """
    #default duration of 15 days
    duration = 15*24*60*60
    #default interval of 6 hrs
    interval = 6*60*60
    #if it is high res
    if(int(ensemble_number) == 52):
        #duration of 10 days
        duration = 10*24*60*60
        #interval of 3 hrs
        #interval = 3*60*60     
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    old_file = open(namelist_file)
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    for line in old_file:
        if line.strip().startswith('BS_opt_Qinit'):
            if (qinit_file):
                new_file.write('BS_opt_Qinit       =.true.\n')
            else:
                new_file.write('BS_opt_Qinit       =.false.\n')
        elif line.strip().startswith('Vlat_file'):
            new_file.write('Vlat_file          =\'../../rapid/input/%s\'\n' % vlat_file)
        elif line.strip().startswith('Qout_file'):
            new_file.write('Qout_file          =\'../../rapid/output/%s\'\n' % qout_file)
        elif line.strip().startswith('ZS_TauM'):
            new_file.write('ZS_TauM            =%s\n' % duration)
        elif line.strip().startswith('ZS_dtM'):
            new_file.write('ZS_dtM             =%s\n' % 86400)
        elif line.strip().startswith('ZS_TauR'):
            new_file.write('ZS_TauR            =%s\n' % interval)
        elif line.strip().startswith('Qinit_file'):
            if (qinit_file):
                new_file.write('Qinit_file          =\'../../rapid/input/%s\'\n' % qinit_file)
            else:
                new_file.write('Qinit_file          =\'\'\n')
        else:
            new_file.write(line)
             
    #close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    #Remove original file
    os.remove(namelist_file)
    #Move new file
    move(abs_path, namelist_file)

def run_RAPID_single_watershed(rapid_files_location, watershed, 
                               forecast_date_timestep, ensemble_number,
                               lock, queue):
    """
    run RAPID on single watershed after ECMWF prepared
    """
    input_directory = os.path.join(rapid_files_location, 'input', 
        watershed)
    #loop through all the rapid_namelist files in directory
    file_list = glob(os.path.join(input_directory,'rapid_namelist_*.dat'))
    for namelist_file in file_list:
        basin_name = os.path.basename(namelist_file)[15:-4]
            
        #make output dir if does not exist
        output_directory = os.path.join(rapid_files_location,
            'output', watershed, forecast_date_timestep)
        try:
            os.makedirs(output_directory)
        except OSError:
            pass

        lock.acquire() #make sure no overlap occurs

        #remove link to old RAPID namelist file
        old_namelist_file = os.path.join(rapid_files_location,
            'run','rapid_namelist')
        try:
            os.unlink(old_namelist_file)
            os.remove(old_namelist_file)
        except OSError:
            pass
        #If any part in here fails, skip and 
        # tell the user what the error is
        try:
            #change the new RAPID namelist file
            print "Updating namelist file for: " + basin_name + " " + str(ensemble_number)
            inflow_file_name = 'm3_riv_bas_'+ str(ensemble_number) + '.nc'
            vlat_fpath = '%s/%s/%s' % (watershed, forecast_date_timestep, inflow_file_name)
                
            qout_fpath = '%s/%s/Qout_%s_%s.nc' % (watershed, forecast_date_timestep, 
                                                  basin_name, ensemble_number)
            #check for qinit file
            past_date = datetime.datetime.strptime(forecast_date_timestep[:11],"%Y%m%d.%H") - \
            datetime.timedelta(0,12*60*60)
            past_hour = '1200' if past_date.hour > 11 else '0'
            past_forecast_date_timestep = "%s.%s" % (past_date.strftime("%Y%m%d"), past_hour)                                     
            qinit_file = os.path.join(input_directory, 'Qinit_file_%s_%s.csv' % (basin_name, past_forecast_date_timestep))                                      
            if os.path.exists(qinit_file):
                qinit_fpath = '%s/Qinit_file_%s_%s.csv' % (watershed, basin_name,
                                                              past_forecast_date_timestep)
            else:
                qinit_fpath = None
                
            #update namelist file for current run    
            update_namelist_file(namelist_file, output_directory,
                    vlat_fpath, qout_fpath, qinit_fpath, ensemble_number)

            #change link to new RAPID namelist file 
            os.symlink(os.path.join(rapid_files_location, 
                'input', watershed, 
                'rapid_namelist_%s.dat' % basin_name).replace("\\","/"),
                os.path.join(rapid_files_location,
                'run', 'rapid_namelist').replace("\\","/")
                )
            
            #run RAPID
            print "Running RAPID for: %s Ensemble: %s" % (basin_name, ensemble_number)

            process = subprocess.Popen([os.path.join(rapid_files_location,'run',
                'run_rapid.sh').replace("\\","/")], shell=True)
                
            sleep(4) #give rapid time to read namelist file

            lock.release()
            queue.put(process)
            process.communicate() #wait for RAPID to finish
        
        except Exception as ex:
            #skip the run
            lock.release()
            print ex
            pass

def prepare_all_inflow_ECMWF(args):
    """
    prepare all ECMWF files for rapid
    """
    lock = args[1]
    queue = args[2]
    watershed = args[0][0]
    forecast = args[0][1]
    rapid_files_location = args[0][2]
    forecast_split = os.path.basename(forecast).split(".")
    forecast_date_timestep = ".".join(forecast_split[:2])
    ensemble_number = int(forecast_split[2])
    directory_path = os.path.join(rapid_files_location, 'input', 
        watershed)

    #create infloe directory if does not exist
    inflow_directory = os.path.join(directory_path, forecast_date_timestep)
    try:
        os.mkdir(inflow_directory)
    except OSError:
        pass
    
    inflow_file_name = 'm3_riv_bas_%s.nc' % ensemble_number
    output_rapid_inflow_file = os.path.join(inflow_directory, 
        inflow_file_name)
    
    #determine weight table from resolution (52 is high resolution)
    weight_table_file = 'weight_low_res.csv'
    if ensemble_number == 52:
        weight_table_file = 'weight_high_res.csv'
    in_weight_table = os.path.join(directory_path, weight_table_file)
    
    time_start_all = datetime.datetime.utcnow()
    try:
        #prepare ECMWF file for RAPID
        #optional argument ... time interval?
        print "Computing Inflow ECMWF for RAPID. Watershed: %s %s %s  ..." % \
            (watershed, forecast_date_timestep, ensemble_number)
        RAPIDinflowECMWF_tool = CreateInflowFileFromECMWFRunoff()
        RAPIDinflowECMWF_tool.execute(forecast, 
            in_weight_table, output_rapid_inflow_file)

        time_end_ecmwf = datetime.datetime.utcnow()

        #run RAPID
        print "Running RAPID. Watershed: %s %s %s  ..." % \
            (watershed, forecast_date_timestep, ensemble_number)
        run_RAPID_single_watershed(rapid_files_location, watershed, 
                               forecast_date_timestep, ensemble_number,
                               lock, queue)
    except Exception as ex:
        print ex
        print "Skipping RAPID analysis for %s %s %s  ..." % \
            (watershed, forecast_date_timestep, ensemble_number)
        return False
    
    #record compute times    
    time_end_all = datetime.datetime.utcnow()
    print "Watershed: %s %s %s Time to Compute: %s ..." % \
        (watershed, forecast_date_timestep, ensemble_number, time_end_all - time_start_all)
    print "Watershed: %s %s %s Time to Downscale ECMWF: %s ..." % \
        (watershed, forecast_date_timestep, ensemble_number, time_end_ecmwf-time_start_all)
    print "Watershed: %s %s %s Time to run RAPID: %s ..." % \
        (watershed, forecast_date_timestep, ensemble_number, time_end_all-time_end_ecmwf)

    return output_rapid_inflow_file
    
            
def make_iterator(args, lock, queue):
    """
    Makes an iterator over args and passes the lock an queue to each element.
    """
    return ((arg, lock, queue) for arg in args)  
    
#------------------------------------------------------------------------------
#main process
#------------------------------------------------------------------------------
#C:\Python27\ArcGIS10.2\python.exe C:\Users\byu_rapid\Documents\RAPID\py_scripts\rapid_process_async.py
if __name__ == "__main__":
    rapid_files_location = '/home/alan/work/rapid'
    ecmwf_forecast_location = "/home/alan/work/ecmwf"
    ckan_api_endpoint = 'http://ciwckan.chpc.utah.edu'
    ckan_api_key = '8dcc1b34-0e09-4ddc-8356-df4a24e5be87'
    download_ecmwf = True
    initialize_flows = True

    #get list of watersheds in rapid directory
    watersheds = [d for d in os.listdir(os.path.join(rapid_files_location,'input')) \
                    if os.path.isdir(os.path.join(rapid_files_location,'input', d))]  
    time_begin_all = datetime.datetime.utcnow()
    date_string = time_begin_all.strftime('%Y%m%d')
    #date_string = datetime.datetime(2015,6,6).strftime('%Y%m%d')

    #download all files for today
    if download_ecmwf:
        ecmwf_folders = ftp_ecmwf_download.download_all_ftp(ecmwf_forecast_location,
                                                            'Runoff.%s*.netcdf.tar.gz' % date_string)
    else:
        ecmwf_folders = glob(os.path.join(ecmwf_forecast_location,
                                          'Runoff.'+date_string+'*.netcdf'))

    #prepare ECMWF files
    time_start_prepare = datetime.datetime.utcnow()
    for ecmwf_folder in ecmwf_folders:
        ecmwf_forecasts = glob(os.path.join(ecmwf_folder,'*.runoff.netcdf'))
        #make the largest files first
        ecmwf_forecasts.sort(key=os.path.getsize, reverse=True)
    
        #loop prepare forecasts async multiple processes at a time if available
        pool = multiprocessing.Pool()
        lock = multiprocessing.Manager().Lock()
        queue = multiprocessing.Manager().Queue()
        combinations = list(itertools.product(watersheds, 
                                              ecmwf_forecasts, 
                                              [rapid_files_location]))
        #chunksize=1 makes it so there is only one task per process
        result = pool.imap(prepare_all_inflow_ECMWF, 
                          make_iterator(combinations, lock, queue),
                          chunksize=1)
        pool.close()
        pool.join()
    
    
        #remove RAPID output files
        rapid_output_files=glob(os.path.join(rapid_files_location,
            'run','*_rapid_stdout.txt'))
        for rapid_output_file in rapid_output_files:
            try:
                os.remove(rapid_output_file)
            except OSError:
                pass
        
        #remove intermediate files
        for ecmwf_folder in ecmwf_folders:
            ecmwf_file = glob(os.path.join(ecmwf_folder,'*.runoff.netcdf'))[0]
            forecast_split = os.path.basename(ecmwf_file).split(".")
            forecast_date_timestep = ".".join(forecast_split[0:2])
            for watershed in watersheds:
                rmtree(os.path.join(rapid_files_location, 'input', 
                    watershed, forecast_date_timestep))

        #initialize flows for next run
        if initialize_flows:
            #create new init flow files
            for watershed in watersheds:
                print "Initializing flows for ", watershed
                input_directory = os.path.join(rapid_files_location, 'input', watershed)
                path_to_watershed_files = os.path.join(rapid_files_location, 'output', watershed)
                forecast_date_timestep = None
                #finds the current output from downscaled ECMWF forecasts
                if os.path.exists(path_to_watershed_files):
                    forecast_date_timestep = sorted([d for d in os.listdir(path_to_watershed_files) \
                                        if os.path.isdir(os.path.join(path_to_watershed_files, d))],
                                         reverse=True)[0]
                if forecast_date_timestep:
                    #loop through all the rapid_namelist files in directory
                    file_list = glob(os.path.join(input_directory,'rapid_namelist_*.dat'))
                    forecast_directory = os.path.join(path_to_watershed_files, forecast_date_timestep)
                    for namelist_file in file_list:
                        basin_name = os.path.basename(namelist_file)[15:-4]
                        basin_files = find_current_rapid_output(forecast_directory, basin_name)
                        compute_initial_rapid_flows(basin_files, basin_name, 
                                                    input_directory, forecast_date_timestep)  
    
    time_finish_prepare = datetime.datetime.utcnow()
 
    #upload new datasets
    data_manager = ECMWFRAPIDDatasetManager(ckan_api_endpoint,
                                            ckan_api_key)
                                            
    data_manager.zip_upload_resources(os.path.join(rapid_files_location, 'output'))
    #delete local datasets
    for item in os.listdir(os.path.join(rapid_files_location, 'output')):
        rmtree(os.path.join(rapid_files_location, 'output', item))

    time_end = datetime.datetime.utcnow()
    #print time to complete all
    print "Time Begin All: ", time_begin_all
    print "Time to Download: ", (time_start_prepare - time_begin_all)
    print "Time Start Prepare: ", time_start_prepare
    print "Time to Prepare: ", (time_finish_prepare - time_start_prepare)
    print "Time Finish All: ", time_end
    print "TOTAL TIME: ", (time_end-time_begin_all)
