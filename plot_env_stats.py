# Load Libraries
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import glob
import os
import copy
import math
import scipy

#Saveadd
#saveaddp = 'all_200'
#saveadd = '_all_200'

saveaddp = 'all_100_p25'
saveadd = '_all_100_p25'

#saveaddp = 'all_100_p25_120'
#saveadd = '_all_100_p25_120'

#saveaddp = 'all_100_p20_120'
#saveadd = '_all_100_p20_120'

#saveaddp = 'all_150'
#saveadd = '_all_150'

#saveaddp = 'all_100_p15_120'
#saveadd = '_all_100_p15_120'

#saveaddp = 'all_100_p10_120'
#saveadd = '_all_100_p10_120'

print(saveadd)
# Create list of cases to loop through
cases = ['ARG1.1-R_old','ARG1.2-R','BRA1.1-R','BRA1.2-R','AUS1.1-R','DRC1.1-R','PHI1.1-R','PHI2.1-R','WPO1.1-R','USA1.1-R','RSA1.1-R']
cases = ['ARG1.1-R','ARG1.2-R','BRA1.1-R','BRA2.1-R','AUS1.1-R','DRC1.1-R','PHI1.1-R','PHI2.1-R','SIO1.1-R','SAU1.1-R','WPO1.1-R','USA1.1-R','USA3.1-R']
#cases = ['WPO1.1-R']
#cases = ['ARG1.1-R_old','ARG1.2-R','BRA1.1-R','BRA1.2-R','AUS1.1-R'] #'DRC1.1-R','PHI1.1-R','PHI1.1-RPR','PHI2.1-R','WPO1.1-R','USA1.1-R','RSA1.1-R','BRA1.1-RPR','DRC1.1-RCN','DRC1.1-RCR']
#cases = ['ARG1.1-R_old']

# Create colors for the different cases
scolors = OrderedDict()
scolors['ARG1.1-R_old'] = 'darkorange'
scolors['ARG1.2-R'] = 'orange'
scolors['PHI2.1-R'] = 'darkviolet'
scolors['BRA1.2-R'] = 'gold'
scolors['AUS1.1-R'] = 'mediumorchid'
scolors['PHI1.1-R'] = 'mediumpurple'
scolors['PHI1.1-RPR'] = 'mediumpurple'
scolors['USA1.1-R'] = 'lightblue'
scolors['BRA1.1-R'] = 'yellowgreen'
scolors['BRA1.1-RPR'] = 'yellowgreen'
scolors['WPO1.1-R'] = 'dodgerblue'
scolors['DRC1.1-R'] = 'green'
scolors['DRC1.1-RCR'] = 'green'
scolors['RSA1.1-R'] = 'olivedrab'

var = OrderedDict()
var[0] = OrderedDict()
var[0]['titlename'] = 'ML CIN'
var[0]['units'] = 'J/kg'
var[0]['varname'] = 'mlcin'
var[0]['varins'] = 'mlcins'
var[0]['varin'] = 'mlcin'
#var[0]['bins'] = np.arange(-300,0.1,2)
var[0]['bins'] = np.arange(0,500,2)
var[0]['sbins'] = np.arange(0,1.01,0.001)

var[1] = OrderedDict()
var[1]['titlename'] = 'ML CAPE'
var[1]['units'] = 'J/kg'
var[1]['varname'] = 'mlcape'
var[1]['varins'] = 'mlcapes'
var[1]['varin'] = 'mlcape'
var[1]['bins'] = np.arange(0,5001,20)
var[1]['sbins'] = np.arange(0,1.01,0.001)

var[2] = OrderedDict()
var[2]['titlename'] = 'Lowlevel Wind Shear'
var[2]['units'] = 'm/s'
var[2]['varname'] = 'Shear_ll'
var[2]['varins'] = 'shr_lls'
var[2]['varin'] = 'shr_ll'
var[2]['bins'] =  np.arange(0,60,0.2)
var[2]['sbins'] = np.arange(0,1.01,0.001)

var[3] = OrderedDict()
var[3]['titlename'] = 'Midlevel Wind Shear'
var[3]['units'] = 'm/s'
var[3]['varname'] = 'Shear_ml'
var[3]['varins'] = 'shr_mls'
var[3]['varin'] = 'shr_ml'
var[3]['bins'] =  np.arange(0,60,0.2)
var[3]['sbins'] = np.arange(0,1.01,0.001)

var[4] = OrderedDict()
var[4]['titlename'] = 'Midlevel Relative Humidity'
var[4]['units'] = 'unitless'
var[4]['varname'] = 'RH_ml'
var[4]['varins'] = 'rh_mls'
var[4]['bins'] = np.arange(0.0,1.01,0.01)
var[4]['sbins'] = np.arange(0,1.01,0.001)

var[5] = OrderedDict()
var[5]['titlename'] = 'Lowlevel Relative Humidity'
var[5]['units'] = 'unitless'
var[5]['varname'] = 'RH_ll'
var[5]['varins'] = 'rh_lls'
var[5]['bins'] = np.arange(0.0,1.01,0.01)
var[5]['sbins'] = np.arange(0,1.01,0.001)

var[6] = OrderedDict()
var[6]['titlename'] = 'Total Column Water Vapor'
var[6]['units'] = 'kg/m^2'
var[6]['varname'] = 'tcwv'
var[6]['varins'] = 'tcwvs'
var[6]['bins'] = np.arange(0,80.1,1)
var[6]['sbins'] = np.arange(0,1.01,0.001)

var[7] = OrderedDict()
var[7]['titlename'] = 'Total Column Condensate'
var[7]['units'] = 'kg/m^2'
var[7]['varname'] = 'tcwc'
var[7]['varins'] = 'tcwcs'
var[7]['bins'] = np.arange(0,100.1,1)
var[7]['sbins'] = np.arange(0,1.01,0.001)

var[8] = OrderedDict()
var[8]['titlename'] = 'Surface RH'
var[8]['units'] = 'unitless'
var[8]['varname'] = 'rhsfc'
var[8]['varins'] = 'rh_sfcs'
var[8]['bins'] = np.arange(0.0,1.11,0.01)
var[8]['sbins'] = np.arange(0,1.01,0.001)

var[9] = OrderedDict()
var[9]['titlename'] = '850hPa RH'
var[9]['units'] = 'unitless'
var[9]['varname'] = 'rh850'
var[9]['varins'] = 'rh_850s'
var[9]['bins'] = np.arange(0.0,1.11,0.01)
var[9]['sbins'] = np.arange(0,1.01,0.001)

var[10] = OrderedDict()
var[10]['titlename'] = '500hPa RH'
var[10]['units'] = 'unitless'
var[10]['varname'] = 'rh500'
var[10]['varins'] = 'rh_500s'
var[10]['bins'] = np.arange(0.0,1.11,0.01)
var[10]['sbins'] = np.arange(0,1.01,0.001)

var[11] = OrderedDict()
var[11]['titlename'] = '250hPa RH'
var[11]['units'] = 'unitless'
var[11]['varname'] = 'rh250'
var[11]['varins'] = 'rh_250s'
var[11]['bins'] = np.arange(0.0,1.11,0.01)
var[11]['sbins'] = np.arange(0,1.01,0.001)

var[12] = OrderedDict()
var[12]['titlename'] = 'Surface Temperature'
var[12]['units'] = 'K'
var[12]['varname'] = 'tsfc'
var[12]['varins'] = 't_sfcs'
var[12]['bins'] = np.arange(0,350.1,0.1)
var[12]['sbins'] = np.arange(0,1.01,0.001)

var[13] = OrderedDict()
var[13]['titlename'] = '850 Temperature'
var[13]['units'] = 'K'
var[13]['varname'] = 't850'
var[13]['varins'] = 't_850s'
var[13]['bins'] = np.arange(0,600.1,0.1)
var[13]['sbins'] = np.arange(0,1.01,0.001)

var[14] = OrderedDict()
var[14]['titlename'] = '500 Temperature'
var[14]['units'] = 'K'
var[14]['varname'] = 't500'
var[14]['varins'] = 't_500s'
var[14]['bins'] = np.arange(0,600.1,0.1)
var[14]['sbins'] = np.arange(0,1.01,0.001)

var[15] = OrderedDict()
var[15]['titlename'] = '250 Temperature'
var[15]['units'] = 'K'
var[15]['varname'] = 't250'
var[15]['varins'] = 't_250s'
var[15]['bins'] = np.arange(0,600.1,0.1)
var[15]['sbins'] = np.arange(0,1.01,0.001)

var[16] = OrderedDict()
var[16]['titlename'] = 'Surface Winds'
var[16]['units'] = 'unitless'
var[16]['varname'] = 'spdsfc'
var[16]['varins'] = 'spd_sfcs'
var[16]['bins'] = np.arange(0.0,50.1,0.5)
var[16]['sbins'] = np.arange(0,1.01,0.001)

var[17] = OrderedDict()
var[17]['titlename'] = '850hPa Winds'
var[17]['units'] = 'unitless'
var[17]['varname'] = 'spd850'
var[17]['varins'] = 'spd_850s'
var[17]['bins'] = np.arange(0.0,50.1,0.5)
var[17]['sbins'] = np.arange(0,1.01,0.001)

var[18] = OrderedDict()
var[18]['titlename'] = '500hPa Winds'
var[18]['units'] = 'unitless'
var[18]['varname'] = 'spd500'
var[18]['varins'] = 'spd_500s'
var[18]['bins'] = np.arange(0.0,50.1,0.5)
var[18]['sbins'] = np.arange(0,1.01,0.001)

var[19] = OrderedDict()
var[19]['titlename'] = '250hPa Winds'
var[19]['units'] = 'unitless'
var[19]['varname'] = 'spd250'
var[19]['varins'] = 'spd_250s'
var[19]['bins'] = np.arange(0.0,50.1,0.5)
var[19]['sbins'] = np.arange(0,1.01,0.001)

var[20] = OrderedDict()
var[20]['titlename'] = 'Surface RV'
var[20]['units'] = 'unitless'
var[20]['varname'] = 'rvsfc'
var[20]['varins'] = 'rv_sfcs'
var[20]['bins'] = np.arange(0.0,0.02,0.001)
var[20]['sbins'] = np.arange(0,1.01,0.001)

var[21] = OrderedDict()
var[21]['titlename'] = '850hPa RV'
var[21]['units'] = 'unitless'
var[21]['varname'] = 'rv850'
var[21]['varins'] = 'rv_850s'
var[21]['bins'] = np.arange(0.0,0.02,0.001)
var[21]['sbins'] = np.arange(0,1.01,0.001)

var[22] = OrderedDict()
var[22]['titlename'] = '500hPa RV'
var[22]['units'] = 'unitless'
var[22]['varname'] = 'rv500'
var[22]['varins'] = 'rv_500s'
var[22]['bins'] = np.arange(0.0,0.02,0.001)
var[22]['sbins'] = np.arange(0,1.01,0.001)

var[23] = OrderedDict()
var[23]['titlename'] = '250hPa RV'
var[23]['units'] = 'unitless'
var[23]['varname'] = 'rv250'
var[23]['varins'] = 'rv_250s'
var[23]['bins'] = np.arange(0.0,0.02,0.001)
var[23]['sbins'] = np.arange(0,1.01,0.001)

var[24] = OrderedDict()
var[24]['titlename'] = 'MidLevel Lapse Rates'
var[24]['units'] = 'K/m'
var[24]['varname'] = 't_mllr'
var[24]['varins'] = 't_mllrs'
var[24]['bins'] = np.arange(2.0,10.,0.1)/1000
var[24]['sbins'] = np.arange(0,1.01,0.001)

var[25] = OrderedDict()
var[25]['titlename'] = 'ML CAPE (GB)'
var[25]['units'] = 'J/kg'
var[25]['varname'] = 'mlcapeGB'
var[25]['varins'] = 'mlcapesGB'
var[25]['varin'] = 'mlcapeGB'
var[25]['bins'] = np.arange(0,5001,20)
var[25]['sbins'] = np.arange(0,1.01,0.001)

var[26] = OrderedDict()
var[26]['titlename'] = 'ML CAPE (GB1)'
var[26]['units'] = 'J/kg'
var[26]['varname'] = 'mlcapeGB1'
var[26]['varins'] = 'mlcapesGB1'
var[26]['varin'] = 'mlcapeGB1'
var[26]['bins'] = np.arange(0,5001,20)
var[26]['sbins'] = np.arange(0,1.01,0.001)

var[27] = OrderedDict()
var[27]['titlename'] = 'MU CAPE (GB)'
var[27]['units'] = 'J/kg'
var[27]['varname'] = 'mucapeGB'
var[27]['varins'] = 'mucapesGB'
var[27]['varin'] = 'mucapeGB'
var[27]['bins'] = np.arange(0,5001,20)
var[27]['sbins'] = np.arange(0,1.01,0.001)

var[28] = OrderedDict()
var[28]['titlename'] = 'ML CIN (GB1)'
var[28]['units'] = 'J/kg'
var[28]['varname'] = 'mlcinGB1'
var[28]['varins'] = 'mlcinsGB1'
var[28]['varin'] = 'mlcinGB1'
var[28]['bins'] = np.arange(0,601,2)
var[28]['sbins'] = np.arange(0,1.01,0.001)

var[29] = OrderedDict()
var[29]['titlename'] = 'ML CIN (GB)'
var[29]['units'] = 'J/kg'
var[29]['varname'] = 'mlcinGB'
var[29]['varins'] = 'mlcinsGB'
var[29]['varin'] = 'mlcinGB'
var[29]['bins'] = np.arange(0,601,2)
var[29]['sbins'] = np.arange(0,1.01,0.001)

var[30] = OrderedDict()
var[30]['titlename'] = 'MU CIN (GB)'
var[30]['units'] = 'J/kg'
var[30]['varname'] = 'mucinGB'
var[30]['varins'] = 'mucinsGB'
var[30]['varin'] = 'mucinGB'
var[30]['bins'] = np.arange(0,601,2)
var[30]['sbins'] = np.arange(0,1.01,0.001)

#subbox = np.array([3,5,9,27])
subbox = np.array([3,4,5,9,27,54])
subbox = np.array([3,4,5,9,27])

pparr = [1,5,10,25,75,90,95,99]
min_val_2d = np.zeros((len(var),len(subbox)))
max_val_2d = np.zeros((len(var),len(subbox)))
p_val_2d = OrderedDict()
for pp in pparr:
    p_val_2d[pp] = np.zeros((len(var),len(subbox)))

#for i in np.arange(0,len(files)):
for v in np.arange(0,len(var)):
    print(v)
    varins = var[v]['varins']
    for s in np.arange(0,len(subbox)):
        var_all_app = []
        for c in np.arange(0,len(cases)):
            cn = cases[c]
            savepath = '/tempest/pmarin//monsoon/ENV/V1/'+cn+'/'+saveaddp+'/'
            files = sorted(glob.glob(savepath+'*SubERA5*.p'))
            for i in np.arange(0,len(files)):
                SEinfilename = savepath+cn+'_Environments_SubERA5'+saveadd+'_'+str(i)+'.p'
                with open(SEinfilename, 'rb') as f:
                    evars = pickle.load(f) # Load environmental variable dictionary
                    var_all = evars[varins]
                    varn = var_all[cn,s]
                var_all_app = np.append(var_all_app,varn)

        min_val_2d[v,s] = np.nanmin(var_all_app)
        max_val_2d[v,s] = np.nanmax(var_all_app)
        for pp in pparr:
            p_val_2d[pp][v,s] = np.nanpercentile(var_all_app,pp)

# Create array for different condensate mixing ratio screens
cscrs = [0.1,0.01,0.001,0.0001,0.00001,0.000001,9999]

# Loop through screens of condensate mixing ratio
for c in np.arange(0,len(cscrs)):
    cscr = cscrs[c]/1000 # kg/kg of condensate
    cstr_txt = 'c'+str(cscrs[c])
    print(cstr_txt)

    #cscr = 9999/1000 # kg/kg of condensate
    #cstr_txt = 'c9999'


    # Create dicrionaty variables for statistics to be saved
    e5_tcwc = OrderedDict()
    avg_save = OrderedDict()
    std_save = OrderedDict()
    ske_save = OrderedDict()
    kur_save = OrderedDict()
    cnt_save = OrderedDict()
    max_save = OrderedDict()
    min_save = OrderedDict()
    p99_save = OrderedDict()
    p95_save = OrderedDict()
    p05_save = OrderedDict()
    p01_save = OrderedDict()
    stdmm_save = OrderedDict()
    stdmma_save = OrderedDict()
    stdp99_save = OrderedDict()
    stdp95_save = OrderedDict()
    stdp90_save = OrderedDict()
    stdp75_save = OrderedDict()
    ran_save = OrderedDict()
    ran0199_save = OrderedDict()
    ran0595_save = OrderedDict()
    ran1090_save = OrderedDict()
    ran2575_save = OrderedDict()
    q1_save = OrderedDict()
    q3_save = OrderedDict()
    qcd_save = OrderedDict()
    cd_save = OrderedDict()
    qcdmma_save = OrderedDict()
    cdmma_save = OrderedDict()
    q1n_save = OrderedDict()
    q3n_save = OrderedDict()
    qcdn_save = OrderedDict()
    cdn_save = OrderedDict()
    q1na_save = OrderedDict()
    q3na_save = OrderedDict()
    qcdna_save = OrderedDict()
    cdna_save = OrderedDict()

    # Loop through cases
    for c in np.arange(0,len(cases)):
    #for c in np.arange(len(cases)-1,len(cases)):
        cn = cases[c]
        print(cn)

        # Data pathnames
        savepath = '/tempest/pmarin/monsoon/ENV/V1/'+cn+'/'+saveaddp+'/'
        if not os.path.exists(savepath):
           os.makedirs(savepath)

        # Location for plotting histograms
        plotpath = '/tempest/pmarin/monsoon/ENV/V1/Plots/'+cn+'/HIST/'
        if not os.path.exists(plotpath):
           os.makedirs(plotpath)

        # Grab all data for each reanalysis box
        files = sorted(glob.glob(savepath+'*SubERA5*.p'))

        # Array of sub-reanalysis box sizes
        subbox = np.array([3,4,5,9,27,54]) # split up 0.25 boxes into 3x3, 4x4, etc...
        subbox = np.array([3,4,5,9,27]) # split up 0.25 boxes into 3x3, 4x4, etc...
        gss_deg = 0.25/subbox
        gss_deg_km = 0.25/subbox*111.000 # calculate ~ size of subreanalysis box sizes

        # Loop through variables
        for v in np.arange(0,len(var)):
            print('Var#:',v)
            titlename = var[v]['titlename']
            units = var[v]['units']
            varname = var[v]['varname']
            varins = var[v]['varins']
            bins = var[v]['bins']
            #print('Analysis Variable:'+varins)

            std_save[cn,v] = np.zeros((len(subbox),len(files)))
            max_save[cn,v] = np.zeros((len(subbox),len(files)))
            min_save[cn,v] = np.zeros((len(subbox),len(files)))
            p99_save[cn,v] = np.zeros((len(subbox),len(files)))
            p95_save[cn,v] = np.zeros((len(subbox),len(files)))
            p05_save[cn,v] = np.zeros((len(subbox),len(files)))
            p01_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdmm_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdmma_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdp99_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdp95_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdp90_save[cn,v] = np.zeros((len(subbox),len(files)))
            stdp75_save[cn,v] = np.zeros((len(subbox),len(files)))
            ran_save[cn,v] = np.zeros((len(subbox),len(files)))
            ran0199_save[cn,v] = np.zeros((len(subbox),len(files)))
            ran0595_save[cn,v] = np.zeros((len(subbox),len(files)))
            ran1090_save[cn,v] = np.zeros((len(subbox),len(files)))
            ran2575_save[cn,v] = np.zeros((len(subbox),len(files)))
            ske_save[cn,v] = np.zeros((len(subbox),len(files)))
            kur_save[cn,v] = np.zeros((len(subbox),len(files)))
            cnt_save[cn,v] = np.zeros((len(subbox),len(files)))
            avg_save[cn,v] = np.zeros((len(subbox),len(files)))
            q1_save[cn,v] = np.zeros((len(subbox),len(files)))
            q3_save[cn,v] = np.zeros((len(subbox),len(files)))
            qcd_save[cn,v] = np.zeros((len(subbox),len(files)))
            cd_save[cn,v] = np.zeros((len(subbox),len(files)))
            qcdmma_save[cn,v] = np.zeros((len(subbox),len(files)))
            cdmma_save[cn,v] = np.zeros((len(subbox),len(files)))

            q1n_save[cn,v] = np.zeros((len(subbox),len(files)))
            q3n_save[cn,v] = np.zeros((len(subbox),len(files)))
            qcdn_save[cn,v] = np.zeros((len(subbox),len(files)))
            cdn_save[cn,v] = np.zeros((len(subbox),len(files)))

            q1na_save[cn,v] = np.zeros((len(subbox),len(files)))
            q3na_save[cn,v] = np.zeros((len(subbox),len(files)))
            qcdna_save[cn,v] = np.zeros((len(subbox),len(files)))
            cdna_save[cn,v] = np.zeros((len(subbox),len(files)))

            # Loop through each reanalysis sized domain for each variable and each case
            for i in np.arange(0,len(files)):
                Einfilename = savepath+cn+'_Environments_ERA5'+saveadd+'.p'
                SEinfilename = savepath+cn+'_Environments_SubERA5'+saveadd+'_'+str(i)+'.p'

                # Load Mean Environmental Variables for ERA5-sized boxes
                with open(Einfilename, 'rb') as f:
                    evar = pickle.load(f) # Load environmental variable dictionary
                    var_E5 = evar[varins]
                    var_tcwc = evar['tcwcs'][cn]
                    var_tcwv = evar['tcwvs'][cn]

                # Load Mean Environmental Variables for Sub ERA5-sized boxes
                with open(SEinfilename, 'rb') as f:
                    evars = pickle.load(f) # Load environmental variable dictionary
                    var_all = evars[varins]
                    # Create condensate screen variables
                    if 'sfc' in varins:
                        sc_scr = evars['c_sfcs']
                    elif '850' in varins:
                        sc_scr = evars['c_850s']
                    elif '500' in varins:
                        sc_scr = evars['c_500s']
                    elif '250' in varins:
                        sc_scr = evars['c_250s']
                    elif '_ml' in varins:
                        sc_scr = evars['c_mls']
                    elif '_ll' in varins:
                        sc_scr = evars['c_lls']
                    elif 'mlcape' in varins:
                        sc_scr = evars['c_lls']
                    elif 'mucape' in varins:
                        sc_scr = evars['c_lls'] 
                    elif 'mlcin' in varins:
                        sc_scr = evars['c_lls']
                    elif 'mucin' in varins:
                        sc_scr = evars['c_lls']
                
                fig,ax = plt.subplots(1,1,figsize=[10,4])
                max_freq = 0
                for s in np.arange(0,len(subbox)):

                    varn = var_all[cn,s]
                    varn = np.abs(varn)

                    # Change data to nan if does not meet screen
                    if ('sfc' in varins) or ('850' in varins) or ('500' in varins) or ('250' in varins) or ('_ll' in varins) or ('_ml' in varins) or ('mlcape' in varins) or ('mucape' in varins) or ('mlcin' in varins) or ('mucin' in varins):
                        scrn = sc_scr[cn,s]
                        #print(np.nanmax(scrn))
                        varn[scrn > cscr] = np.nan

                    #print(np.nanmean(varn))
                    avg_save[cn,v][s,i] = np.nanmean(varn)
                    std_save[cn,v][s,i] = np.nanstd(varn)
                    ran_save[cn,v][s,i] = np.nanmax(varn)-np.nanmin(varn)
                    max_save[cn,v][s,i] = np.nanmax(varn)
                    p99_save[cn,v][s,i] = np.nanpercentile(varn,99)
                    p01_save[cn,v][s,i] = np.nanpercentile(varn,1)
                    p95_save[cn,v][s,i] = np.nanpercentile(varn,95)
                    p05_save[cn,v][s,i] = np.nanpercentile(varn,5)
                    ran0199_save[cn,v][s,i] = p99_save[cn,v][s,i] - p01_save[cn,v][s,i]
                    ran0595_save[cn,v][s,i] = p95_save[cn,v][s,i] - p05_save[cn,v][s,i]
                    ran1090_save[cn,v][s,i] = np.nanpercentile(varn,90) - np.nanpercentile(varn,10)
                    ran2575_save[cn,v][s,i] = np.nanpercentile(varn,75) - np.nanpercentile(varn,25)
                    min_save[cn,v][s,i] = np.nanmin(varn)
                    ske_save[cn,v][s,i] = scipy.stats.skew(varn)
                    kur_save[cn,v][s,i] = scipy.stats.kurtosis(varn)
                    cnt_save[cn,v][s,i] = np.count_nonzero(~np.isnan(varn))
                    q1_save[cn,v][s,i] = np.nanpercentile(varn,25)
                    q3_save[cn,v][s,i] = np.nanpercentile(varn,75)
                    qcd_save[cn,v][s,i] = (q3_save[cn,v][s,i]-q1_save[cn,v][s,i])/(q3_save[cn,v][s,i]+q1_save[cn,v][s,i])
                    cd_save[cn,v][s,i] = std_save[cn,v][s,i]/avg_save[cn,v][s,i]

                    # Local Simulation Min-max scaled statistics
                    varn_mm = (varn-np.nanmin(varn))/(np.nanmax(varn)-np.nanmin(varn))
                    stdmm_save[cn,v][s,i] = np.nanstd(varn_mm)

                    # Global Min-max scaled statistics
                    varn_mma = (varn-min_val_2d[v,s])/(max_val_2d[v,s]-min_val_2d[v,s])
                    stdmma_save[cn,v][s,i] = np.nanstd(varn_mma)
                    cdmma_save[cn,v][s,i] = np.nanstd(varn_mma)/np.nanmean(varn_mma)
                    qcdmma_save[cn,v][s,i] = (np.nanpercentile(varn_mma,75)-np.nanpercentile(varn_mma,25))/(np.nanpercentile(varn_mma,75)+np.nanpercentile(varn_mma,25))

                    varn_p99 = (varn-p_val_2d[1][v,s])/(p_val_2d[99][v,s]-p_val_2d[1][v,s])
                    stdp99_save[cn,v][s,i] = np.nanstd(varn_p99)

                    varn_p95 = (varn-p_val_2d[5][v,s])/(p_val_2d[95][v,s]-p_val_2d[5][v,s])
                    stdp95_save[cn,v][s,i] = np.nanstd(varn_p95)

                    varn_p90 = (varn-p_val_2d[10][v,s])/(p_val_2d[90][v,s]-p_val_2d[10][v,s])
                    stdp90_save[cn,v][s,i] = np.nanstd(varn_p90)

                    varn_p75 = (varn-p_val_2d[25][v,s])/(p_val_2d[75][v,s]-p_val_2d[25][v,s])
                    stdp75_save[cn,v][s,i] = np.nanstd(varn_p75)


                    # Local Simulation Min-scaled statistics
                    q1n_save[cn,v][s,i] = np.nanpercentile(varn-np.nanmin(varn),25)
                    q3n_save[cn,v][s,i] = np.nanpercentile(varn-np.nanmin(varn),75)
                    qcdn_save[cn,v][s,i] = (q3n_save[cn,v][s,i]-q1n_save[cn,v][s,i])/(q3n_save[cn,v][s,i]+q1n_save[cn,v][s,i])
                    cdn_save[cn,v][s,i] = np.nanstd(varn-np.nanmin(varn))/np.nanmean(varn-np.nanmin(varn))

                    # Global Min-scaled statistics
                    q1na_save[cn,v][s,i] = np.nanpercentile(varn-min_val_2d[v,s],25)
                    q3na_save[cn,v][s,i] = np.nanpercentile(varn-min_val_2d[v,s],75)
                    qcdna_save[cn,v][s,i] = (q3na_save[cn,v][s,i]-q1na_save[cn,v][s,i])/(q3na_save[cn,v][s,i]+q1na_save[cn,v][s,i])
                    cdna_save[cn,v][s,i] = np.nanstd(varn-min_val_2d[v,s])/np.nanmean(varn-min_val_2d[v,s])

                    freq,bins = np.histogram(varn,bins=bins)

                    freq = np.insert(freq,0,0)
                    ax.step(bins,freq/np.nansum(freq),label=str(np.round(gss_deg_km[s],1))+' ('+str(np.round(std_save[cn,v][s,i],2))+','+str(int(np.sum(freq)))+')')
                    max_temp = np.nanmax(freq/np.nansum(freq))
                    if max_temp > max_freq:
                        max_freq = copy.deepcopy(max_temp)

                print(len(var_E5[cn]),i)
                ax.plot([var_E5[cn][i],var_E5[cn][i]],[0.0,1.0],'-k')
                ax.set_ylim([0,max_freq+0.01])
                ax.set_title(cn+': Grid'+str(i)+' '+varname)
                ax.set_xlabel(titlename+' ('+units+')')
                ax.set_ylabel('Relative Frequency')
                ax.legend()
                ax.grid()
                plt.tight_layout()
                plt.savefig(plotpath+cn+'_'+varname+str(i)+'_'+cstr_txt+saveadd+'.png')
                plt.close(fig)


        #e5_tcwc[cn,v] = np.reshape(np.tile(np.log(var_tcwc),np.shape(subbox)),np.shape(avg_save[cn,v]))
        #fig,ax = plt.subplots(1,1,figsize=[10,4])
        #a = ax.scatter(np.transpose(np.tile(gss_deg_km,[np.shape(ran_save[cn,v])[1],1])),ran_save[cn,v],c=e5_tcwc[cn,v])
        #cbar = plt.colorbar(a,ax=ax)
        #cbar.ax.set_ylabel('Total Column Condensate (log e)')
        #ax.plot(gss_deg_km,np.nanmean(ran_save[cn,v],axis=1),'-k')
        #ax.set_xscale('log')
        #ax.set_xticks([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xlabel('dx,dy (km)')
        #ax.set_ylabel('Range'+' ('+units+')')
        #ax.set_title(cn+': '+titlename)
        #ax.grid()
        #plt.tight_layout()
        #plt.savefig(plotpath+varname+'_Range_v_Spacingkm_'+cstr_txt+saveadd+'.png')
        #plt.savefig(plotpath+varname+'_Range_v_Spacingkm_'+cstr_txt+saveadd+'.pdf')
        #plt.close(fig)

        #fig,ax = plt.subplots(1,1,figsize=[10,4])
        #a = ax.scatter(np.transpose(np.tile(gss_deg_km,[np.shape(std_save[cn,v])[1],1])),std_save[cn,v],c=e5_tcwc[cn,v])
        #cbar = plt.colorbar(a,ax=ax)
        #cbar.ax.set_ylabel('Total Column Condensate (log e)')
        #ax.plot(gss_deg_km,np.nanmean(std_save[cn,v],axis=1),'-k')
        #ax.set_xscale('log')
        #ax.set_xticks([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xlabel('dx,dy (km)')
        #ax.set_ylabel('Standard Deviation'+' ('+units+')')
        #ax.set_title(cn+': '+titlename)
        #ax.grid()
        #plt.tight_layout()
        #plt.savefig(plotpath+varname+'_STD_v_Spacingkm_'+cstr_txt+saveadd+'.png')
        #plt.savefig(plotpath+varname+'_STD_v_Spacingkm_'+cstr_txt+saveadd+'.pdf')
        #plt.close(fig)

        #fig,ax = plt.subplots(1,1,figsize=[10,4])
        #a = ax.scatter(np.transpose(np.tile(gss_deg_km,[np.shape(avg_save[cn,v])[1],1])),avg_save[cn,v],c=e5_tcwc[cn,v])
        #cbar = plt.colorbar(a,ax=ax)
        #cbar.ax.set_ylabel('Total Column Condensate (log e)')
        #ax.plot(gss_deg_km,np.nanmean(avg_save[cn,v],axis=1),'-k')
        #ax.set_xscale('log')
        #ax.set_xticks([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10,15])
        #ax.set_xlabel('dx,dy (km)')
        #ax.set_ylabel('Mean'+' ('+units+')')
        #ax.set_title(cn+': '+titlename)
        #ax.grid()
        #plt.tight_layout()
        #plt.savefig(plotpath+varname+'_AVG_v_Spacingkm_'+cstr_txt+saveadd+'.png')
        #plt.savefig(plotpath+varname+'_AVG_v_Spacingkm_'+cstr_txt+saveadd+'.pdf')
        #plt.close(fig)

    # open a file, where you ant to store the data
    file = open('/tempest/pmarin/monsoon/ENV/V1/Var_Stats_'+cstr_txt+saveadd+'.p', 'wb')

    data = [avg_save, std_save, ran_save, ran0199_save, ran0595_save, ran1090_save, ran2575_save, ske_save, kur_save, cnt_save, min_save, p01_save, p05_save, p95_save, p99_save, max_save, q1_save, q3_save, qcd_save, cd_save, stdmm_save, stdmma_save, stdp99_save, stdp95_save, stdp90_save, stdp75_save, cdmma_save, qcdmma_save]
    # dump information to that file
    pickle.dump(data, file)

    # close the filee
    file.close()

