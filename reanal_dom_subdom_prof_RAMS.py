# Import Python Libraries
import numpy as np
from datetime import datetime
import sys
import glob
#import rams_tools
import os
import h5py
from collections import OrderedDict
import hdf5plugin
import copy
import pickle
import matplotlib.pyplot as plt
import matplotlib

# Function to decompose domains into reanalysis-sized areas
def grid_decomp(filein,ex_pts,grid_spacing):
    # Grab initial file at the beginning of the analysis period
    rams_file = h5py.File(filein, 'r')

    # Grab lat lon variables from case initial file
    lat = np.array(rams_file['GLAT'])
    lon = np.array(rams_file['GLON'])
    dom_rat = np.shape(lon)[0]/np.shape(lon)[1]

    # Calculate domain size
    ny = np.shape(lat)[0]
    nx = np.shape(lat)[1]

    # Calculate lat lon bounds
    lat_avg = np.nanmean(lat,axis=1)
    lon_avg = np.nanmean(lon,axis=0)

    # Ignore 100 points near grid boundaries ((exclude 10 km near the boundaries))
    lat0 = lat_avg[ex_pts]
    lat1 = lat_avg[ny-ex_pts]
    lon0 = lon_avg[ex_pts]
    lon1 = lon_avg[nx-ex_pts]

    lon_arr = np.arange(lon0,lon1+0.01,grid_spacing)
    lat_arr =  np.arange(lat0,lat1+0.01,grid_spacing)

    return lon_avg,lat_avg,lon_arr,lat_arr,lon0,lat0,lon1,lat1,nx,ny


# Loop through cases simulations
cases = ['SIO1.1-R','ARG1.1-R','ARG1.2-R','PHI2.1-R']
cases = ['BRA2.1-R','SAU1.1-R','USA3.1-R']
savepathd = '/tempest/pmarin/monsoon/ENV/V1/Doms/' # Savepath for saving domain images
savepathp = '/tempest/pmarin/monsoon/ENV/V1/' # Savepath for saving profile fidata



# Specify domain features
gs_reanal = 0.25 # Degrees for reanalysis grid boxes
ex_pts = 100 # exlcude # grid points near edges
savename = 'all_100_p25_120'
fid = 119

savename = 'all_100_p25'
fid = 0

#ex_pts = 200 # exlcude # grid points near edges
#savename = 'all_200'

# Define variables for 27KM_AREAS and 1KM_SUBGRID
prof = OrderedDict()
latlon = OrderedDict()
numg = OrderedDict()
profs = OrderedDict()
latlons_id = OrderedDict()
latlons = OrderedDict()
numgs = OrderedDict()

# Create maps of ERA5 decomposition for saving
# Loop through cases
plt.rcParams.update({'font.size': 16})
for c in np.arange(0,len(cases)):
    cn = cases[c]
    print(cn)
    continue
    prof[cn] = OrderedDict()
    profs[cn] = OrderedDict()
    # get all lite files
    files = sorted(glob.glob('/monsoon/MODEL/LES_MODEL_DATA/V1/'+cn+'-V1/G3/out_30s/a-L*g3.h5'))

    filein = files[fid] # Choose first lite file

    # Calculate lat /lon bounds
    lon_avg,lat_avg,lon_arr,lat_arr,lon0,lat0,lon1,lat1,nx,ny = grid_decomp(filein,ex_pts,gs_reanal)

    # Variable for plotting
    rams_file = h5py.File(filein, 'r')

    # Grab lat lon variables from case initial file
    lat = np.array(rams_file['GLAT'])
    lon = np.array(rams_file['GLON'])
    dom_rat = np.shape(lon)[0]/np.shape(lon)[1]
    
    #plt_var = rams_file['RV'][51,:,:]*1000 # Plot variable

    plt_var = np.sqrt(np.power(rams_file['UP'][1,:,:],2) + np.power(rams_file['VP'][1,:,:],2)) # Plot variable
    
    rams_file.close() # Close file
    
    w_lvls = 100 # Specify number of contour levels for plotting
    # Make Plot
    fig,ax = plt.subplots(1,1,figsize=[7.5,7*dom_rat])
    a = ax.contourf(lon,lat,plt_var,levels=w_lvls,cmap=plt.cm.viridis,extend='both')
    for lo in np.arange(0,len(lon_arr)):
        plt.plot([lon_arr[lo],lon_arr[lo]],[np.min(lat_arr),np.max(lat_arr)],'-w')

    for la in np.arange(0,len(lat_arr)):
        plt.plot([np.min(lon_arr),np.max(lon_arr)],[lat_arr[la],lat_arr[la]],'-w')

    # Number reanalysis sized boxes
    cntp = 0    
    # loop through lat and lon regions (add labels)
    for lo in np.arange(0,len(lon_arr)-1):
        for la in np.arange(0,len(lat_arr)-1):

            lon_id0 = np.where(np.abs(lon_avg-lon_arr[lo]) == np.min(np.abs(lon_avg-lon_arr[lo])))[0][0]
            lon_id1 = np.where(np.abs(lon_avg-lon_arr[lo+1]) == np.min(np.abs(lon_avg-lon_arr[lo+1])))[0][0]
            lat_id0 = np.where(np.abs(lat_avg-lat_arr[la]) == np.min(np.abs(lat_avg-lat_arr[la])))[0][0]
            lat_id1 = np.where(np.abs(lat_avg-lat_arr[la+1]) == np.min(np.abs(lat_avg-lat_arr[la+1])))[0][0]
            lon_plot = np.nanmean([lon_avg[lon_id0],lon_avg[lon_id1]])
            lat_plot = np.nanmean([lat_avg[lat_id0],lat_avg[lat_id1]])
            plt.text(lon_plot-0.05,lat_plot-0.05,str(cntp),color='w',fontsize=10)
            cntp = cntp + 1

    plt.grid()
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cbar = plt.colorbar(a,ax=ax)

    #cbar.ax.set_ylabel('Column Maximum Vertical Velocity (m/s)')
    #cbar.ax.set_ylabel('5kmAGL Vapor Mixing Ratio (g/kg)')
    #cbar.ax.set_ylabel('Near Surface Potential Temperature (K)')
    cbar.ax.set_ylabel('Near Surface Wind Speed (m/s)')

    savefile = savepathd+'G3_5KMRV_Domain_Decomp_0p25deg_'+cn+savename+'.pdf'
    fig.savefig(savefile)
    savefile = savepathd+'G3_5KMRV_Domain_Decomp_0p25deg_'+cn+savename+'.png'
    fig.savefig(savefile)

# Create profile variables for saving 
#27KM_AREA Profiles
prof = OrderedDict() # Variables
latlon = OrderedDict() # Lat/lon information
numg = OrderedDict() # Number of Grids
#1KM_SUBGRID_PROFILES
profs = OrderedDict()
latlons_id = OrderedDict()
latlons = OrderedDict()
numgs = OrderedDict()

#Loop through simulation cases
beg_time2 = datetime.now()
for c in np.arange(0,len(cases)):
    cn = cases[c]

    print(cn)
    prof[cn] = OrderedDict()
    profs[cn] = OrderedDict()
    # get all lit files
    files = sorted(glob.glob('/monsoon/MODEL/LES_MODEL_DATA/V1/'+cn+'-V1/G3/out_30s/a-L*g3.h5'))
    # Grab initial file at the beginning of the analysis period   
    filein = files[fid] # Choose first lite file
    #filein = files[120] # Choose first lite file
    print(filein)
    
    # Calculate lat /lon bounds
    lon_avg,lat_avg,lon_arr,lat_arr,lon0,lat0,lon1,lat1,nx,ny = grid_decomp(filein,ex_pts,gs_reanal)
    
    # Grab initial file at the beginning of the analysis period
    rams_file = h5py.File(filein, 'r')
    lat = rams_file['GLAT'][:]
    lon = rams_file['GLON'][:]

    # Variables to pull and grab profiles of
    varns = ['topt','UP','VP','RV','WP','THETA','temp','pres','RTC','SPD','RTL','RTI','RTRI']
    
    #loop through variables
    for v in np.arange(0,len(varns)):
        var_time1 = datetime.now()
        #print(varns[v])
        prof[cn][varns[v]] = OrderedDict()
        profs[cn][varns[v]] = OrderedDict()
        cntp = 0    
        # loop through lat and lon regions
        for lo in np.arange(0,len(lon_arr)-1):
            for la in np.arange(0,len(lat_arr)-1):
                
                # Find ids for each ERA5 subgrid
                lon_id0 = np.where(np.abs(lon_avg-lon_arr[lo]) == np.min(np.abs(lon_avg-lon_arr[lo])))[0][0]
                lon_id1 = np.where(np.abs(lon_avg-lon_arr[lo+1]) == np.min(np.abs(lon_avg-lon_arr[lo+1])))[0][0]
                lat_id0 = np.where(np.abs(lat_avg-lat_arr[la]) == np.min(np.abs(lat_avg-lat_arr[la])))[0][0]
                lat_id1 = np.where(np.abs(lat_avg-lat_arr[la+1]) == np.min(np.abs(lat_avg-lat_arr[la+1])))[0][0]

                #print(len(np.arange(lon_id0,lon_id1)))
                #print(len(np.arange(lat_id0,lat_id1)))
                #print(lon_id0,lon_id1)
                #print(lat_id0,lat_id1)
                #print(cntp,varns[v],filein)

                #Some variables need additional calculations from raw model data
                if varns[v] == 'temp':
                    cp = 1004; # J/kg/K
                    temp_arr = rams_file['THETA'][:,lat_id0:lat_id1,lon_id0:lon_id1] * (rams_file['PI'][:,lat_id0:lat_id1,lon_id0:lon_id1] / cp)
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'pres':
                    cp = 1004; # J/kg/K
                    rd = 287; # J/kg/K
                    p00 = 100000; # Reference Pressure
                    temp_arr = p00 * np.power( (rams_file['PI'][:,lat_id0:lat_id1,lon_id0:lon_id1] / cp), cp/rd)
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'RTC':
                    temp_arr = rams_file['RTP'][:,lat_id0:lat_id1,lon_id0:lon_id1] - rams_file['RV'][:,lat_id0:lat_id1,lon_id0:lon_id1]
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'RTL':
                    temp_arr = rams_file['RCP'][:,lat_id0:lat_id1,lon_id0:lon_id1] + rams_file['RRP'][:,lat_id0:lat_id1,lon_id0:lon_id1] + rams_file['RDP'][:,lat_id0:lat_id1,lon_id0:lon_id1] 
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'RTI':
                    temp_arr = rams_file['RPP'][:,lat_id0:lat_id1,lon_id0:lon_id1] + rams_file['RAP'][:,lat_id0:lat_id1,lon_id0:lon_id1] + rams_file['RSP'][:,lat_id0:lat_id1,lon_id0:lon_id1] 
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'RTRI':
                    temp_arr = rams_file['RHP'][:,lat_id0:lat_id1,lon_id0:lon_id1] + rams_file['RGP'][:,lat_id0:lat_id1,lon_id0:lon_id1] 
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'SPD':
                    temp_arr = np.sqrt(np.power(rams_file['UP'][:,lat_id0:lat_id1,lon_id0:lon_id1],2.0) + np.power(rams_file['VP'][:,lat_id0:lat_id1,lon_id0:lon_id1],2.0))
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))
                elif varns[v] == 'topt':
                    temp_arr = rams_file['TOPT'][lat_id0:lat_id1,lon_id0:lon_id1]
                    avg_prof = np.nanmean(temp_arr)
                else:
                    temp_arr = rams_file[varns[v]][:,lat_id0:lat_id1,lon_id0:lon_id1]
                    avg_prof = np.nanmean(temp_arr,axis=(1,2))

                # Save mean profile for variable within ERA5-sized domains (27KM_AREAS)
                prof[cn][varns[v]][cntp] = copy.deepcopy(avg_prof)

                # Create storage array for sub-ERA5 domains (1KM_SUBGRID)
                profs[cn][varns[v]][cntp] = OrderedDict()
                # Only need to grab the lat lon bounds for one of the variables
                if v == 0:
                    latlon[cn,cntp] = copy.deepcopy(np.array([lon_id0,lon_id1,lat_id0,lat_id1]))
                
                # Sub Reanalysis-domains (i.e., 1KM_SUBGRID)
                sub_era5 = [3,4,5,9,27] # 9km, 7km, 5km, 3km, 1km # Tested various XKM_SUBGRIDS
                # loop through different subgrid sizes, and perform similar domain decomposition and calculating variables
                start_subtime = datetime.now()
                for s in np.arange(0,len(sub_era5)):
                    
                    profs[cn][varns[v]][cntp][sub_era5[s]] = OrderedDict()

                    # Calculate lat lon bounds for subdomains
                    lat_avgs = np.nanmean(lat[lat_id0:lat_id1,lon_id0:lon_id1],axis=1)
                    lon_avgs = np.nanmean(lon[lat_id0:lat_id1,lon_id0:lon_id1],axis=0)
                    #print(lat_avgs)
                    #print(lon_avgs)

                    # Get lat lon bounds from new subdomain
                    lat0s = np.nanmin(lat_avgs)
                    lat1s = np.nanmax(lat_avgs)
                    lon0s = np.nanmin(lon_avgs)
                    lon1s = np.nanmax(lon_avgs)                   
                    
                    # Calculate lat lon bound arrays
                    dd = sub_era5[s]
                    lon_arr_s = np.arange(lon0s,lon1s+0.001,(lon1s-lon0s)/dd)
                    lat_arr_s = np.arange(lat0s,lat1s+0.001,(lat1s-lat0s)/dd)
                    
                    #print(len(lon_arr_s)-1,len(lat_arr_s)-1)
                    
                    cntps = 0                    
                    # loop through lat and lon regions
                    for los in np.arange(0,len(lon_arr_s)-1):
                        for las in np.arange(0,len(lat_arr_s)-1):
                            subloop_time1 = datetime.now()

                            # Find ids for each ERA5 subgrid
                            lon_id0s = np.where(np.abs(lon_avgs-lon_arr_s[los]) == np.min(np.abs(lon_avgs-lon_arr_s[los])))[0][0]
                            lon_id1s = np.where(np.abs(lon_avgs-lon_arr_s[los+1]) == np.min(np.abs(lon_avgs-lon_arr_s[los+1])))[0][0]
                            lat_id0s = np.where(np.abs(lat_avgs-lat_arr_s[las]) == np.min(np.abs(lat_avgs-lat_arr_s[las])))[0][0]
                            lat_id1s = np.where(np.abs(lat_avgs-lat_arr_s[las+1]) == np.min(np.abs(lat_avgs-lat_arr_s[las+1])))[0][0]

                            #print('where',datetime.now()-subloop_time1)
                            #print(len(np.arange(lon_id0s,lon_id1s)))
                            #print(len(np.arange(lat_id0s,lat_id1s)))
                            #print(lon_id0s,lon_id1s)
                            #print(lat_id0s,lat_id1s)
                           
                            lat_idn = np.arange(lat_id0+lat_id0s,lat_id0+lat_id1s)
                            lon_idn = np.arange(lon_id0+lon_id0s,lon_id0+lon_id1s)
                            #print('dims',datetime.now()-subloop_time1)
                            #print(cntps,varns[v],filein)

                            # Some variables need additional calculations
                            if varns[v] == 'temp':
                                cp = 1004; # J/kg/K
                                temp_arr = rams_file['THETA'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] * (rams_file['PI'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] / cp)
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'pres':
                                cp = 1004; # J/kg/K
                                rd = 287; # J/kg/K
                                p00 = 100000; # Reference Pressure
                                temp_arr = p00 * np.power( (rams_file['PI'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] / cp), cp/rd)
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'RTC':
                                temp_arr = rams_file['RTP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] - rams_file['RV'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s]
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'RTL':
                                temp_arr = rams_file['RCP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] + rams_file['RRP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] + rams_file['RDP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s]
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'RTI':
                                temp_arr = rams_file['RPP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] + rams_file['RAP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] + rams_file['RSP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s]
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'RTRI':
                                temp_arr = rams_file['RHP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] + rams_file['RGP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s] 
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'SPD':
                                temp_arr = np.sqrt(np.power(rams_file['UP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s],2.0) + np.power(rams_file['VP'][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s],2.0))
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            elif varns[v] == 'topt':
                                temp_arr = rams_file['TOPT'][lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s]
                                avg_prof = np.nanmean(temp_arr)
                            else:
                                temp_arr = rams_file[varns[v]][:,lat_id0+lat_id0s:lat_id0+lat_id1s,lon_id0+lon_id0s:lon_id0+lon_id1s]
                                avg_prof = np.nanmean(temp_arr,axis=(1,2))
                            
                            #print('read+avg',datetime.now()-subloop_time1)
#                            profs[cntp,cntps] = copy.deepcopy(var_prof)
#                            profs[cn][varns[v]][sub_era5[s]][cntp,cntps] =  
                            profs[cn][varns[v]][cntp][sub_era5[s]][cntps] = copy.deepcopy(avg_prof)                            
                            #print('copy',datetime.now()-subloop_time1)
                            #print('profs',cn,dd,varns[v],cntp,cntps)
                            # Only need to grab the lat lon bounds for one of the variables
                            if v == 0:
                                latlons_id[cn,dd,cntp,cntps] = copy.deepcopy(np.array([np.min(lon_idn),np.max(lon_idn),np.min(lat_idn),np.max(lat_idn)]))
                                latlons[cn,dd,cntp,cntps] = copy.deepcopy(np.array([lon_avg[np.min(lon_idn)],lon_avg[np.max(lon_idn)],lat_avg[np.min(lat_idn)],lat_avg[np.max(lat_idn)]]))
                            cntps = cntps + 1
                            #print('Time to compute 1 subloop:',datetime.now()-subloop_time1)
                    #print(s,dd,cntps)
                    numgs[cn,cntp,s] = cntps
                print(cntp)
                cntp = cntp + 1
        numg[cn] = cntp
 
    rams_file.close()     
    
    # Save variables for each simulation/case 
    savefile = savepathp+'Sim_Profs_range_'+savename+'_'+cn+'.p'
    with open(savefile, 'wb') as file:
        # A new file will be created
        pickle.dump([prof,latlon,numg,profs,latlons,latlons_id,numgs], file)

