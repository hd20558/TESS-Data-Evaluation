# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 12:11:51 2022

@author: Louis Eddershaw
"""
import lightkurve as lk
import lightkurve.periodogram
import matplotlib.pyplot as plt
import numpy as np

figure_size = (15,5)
BINNING = 0.05
detailed_analysis_threshold = 34
y_limits = False #[0.95,1.05] #False
                #       1               2           3               4           5               6           7               8           9           10          11              12              13          14              15              16              17          18          19              20              21          22              23          24          25          26          27              28          29              30          31              32          33
observations = [[274039311,31],[126947245,1],[341849173,7],[441461124,2],[273985862,1],[266980320,1],[301407485,31],[9967126,4],[120896927,4],[231702397,1],[279614617,9],[103633672,14],[413248763,8],[441791294,14],[20291794,23],[321103619,20],[110996418,10],[149950920,21],[55525572,4],[207468071,23],[183596242,3],[350618622,2],[144539611,4],[11465798,9],[1030783,4],[7809321,4],[20448500,9],[23609565,9],[14336130,20],[441420236,1],[136916387,12],[169226822,9],[160708862,38]]
desired_indexes = [index for index in range(len(observations))]
#desired_indexes = [-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]
flatten_indexes = []

observations_to_analyse = [observations[desired_index] for desired_index in desired_indexes]
#print(observations_to_analyse)

for i in range(len(observations_to_analyse)):
    target_identifier = 'TIC {0}'.format(observations_to_analyse[i][0])
    sector_number = observations_to_analyse[i][1]
    print(i + 1, 'out of', len(observations_to_analyse))

    search_result = lk.search_lightcurve(target_identifier, mission='TESS', sector=sector_number, 
                                         exptime=120, author='SPOC')
    
    lc = search_result.download()
    print(target_identifier, 'sector', sector_number, 'lightcurve downloaded')
    #print(lc.keys()) # This will tell you what objects are in the light curve file we have just downloaded
    #print('---------')
    #print(lc.meta.keys()) # This will tell you other keys that are contained in the metadat of the file
    
    # We now need to normalize our light curve so we can plot it and examine the target
    lc = lc.normalize()
    print(target_identifier, 'sector', sector_number, 'lightcurve normalised')
    
    lc = lc.remove_outliers(sigma=10) # clip the light curve to remove deviant datapoints greater than 10 sigma
    mjd_time = lc.time.mjd - lc.time[0].mjd
    # You can play around with the sigma value - what does this do to the light curve?
    # This also removes nans
    #print(lc.flux.shape, lc.flux[0:10])
    print(target_identifier, 'sector', sector_number, 'lightcurve outliers removed')
    
    fig = plt.figure(figsize=figure_size)
    plt.errorbar(mjd_time, lc.flux, lc.flux_err.data, fmt='.', ecolor='LightGrey',label='Unbinned Data')
    # Note the values are called from the associated lightcurve result we defined above
    # x = the time array (in units of MJD) minus the time of the first data point
    # y = flux
    # err = flux_err
    #print(target_identifier, 'sector', sector_number, 'lightcurve unbinned data plotted')
    
    # Looking through the documentation for Lightkurve we can see it has a binning function
    # https://docs.lightkurve.org/reference/api/lightkurve.LightCurve.bin.html#lightkurve.LightCurve.bin
    # Let's bin up the light curve and overplot that on the full data
    if len(observations_to_analyse) < detailed_analysis_threshold or BINNING > 0.05:
        lc_bin = lc.bin(time_bin_size = BINNING) # Default time is in days
        mjd_time_binned = lc_bin.time.mjd - lc_bin.time[0].mjd
        print(target_identifier, 'sector', sector_number, 'lightcurve data binned')
        plt.errorbar(mjd_time_binned, lc_bin.flux, lc_bin.flux_err.data, fmt='.', ecolor='green',label='Binned Data')
        #print(target_identifier, 'sector', sector_number, 'lightcurve binned data plotted')
        plt.title('{0}, Sector {1}, Binning {2} days'.format(target_identifier,sector_number,BINNING))
        plt.legend()
    else:
        plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
    plt.xlabel('Time [days]')
    plt.ylabel('Normalized Flux')
    if y_limits != False:
        plt.ylim(y_limits[0],y_limits[1])
    plt.show()
    
    # We can calculate the standard deviation of the data and compare it to the uncertainty on the flux given
    scatter = np.std(lc.flux) 
    print('Standard deviation = ', scatter)
    print('Mean err = ', np.mean(lc.flux_err))
    
    
    fig = plt.figure(figsize=figure_size)
    pg = lk.periodogram.BoxLeastSquaresPeriodogram.from_lightcurve(lc)
    print(target_identifier, 'sector', sector_number, 'periodogram generated')
    #pg_bin = lk.periodogram.BoxLeastSquaresPeriodogram.from_lightcurve(lc_bin)
    plt.plot(pg.period, pg.power, label='Unbinned Data')
    #plt.plot(pg_bin.period, pg_bin.power, label='Binned Data')
    plt.xlabel('Period [Days]')
    plt.ylabel('BLS Power')
    plt.title('{0}, Sector {1} Periodogram'.format(target_identifier,sector_number))
    plt.legend()
    plt.show()
    
       
    #print(lc.time.mjd.shape, lc.flux.shape, lc.flux_err.shape) #look at how many points are in your lightcurve
    #cut = 2000 # Select a range or part of the lightcurve to look at
    #fig = plt.figure(figsize=figure_size)
    #plt.plot(lc.time.mjd[0:cut] - lc.time[0].mjd,lc.flux[0:cut])
    #print(target_identifier, 'sector', sector_number, 'lightcurve plotted')
    
    #plt.xlabel('Time [days]')
    #plt.ylabel('Normalised Flux')
    #plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
    #if y_limits != False:
    #    plt.ylim(y_limits[0],y_limits[1])
    #plt.show()
    
    
    
    if i in flatten_indexes:
        # Let's define a new lightcurve ojbect where we have flattened the light curve
        lc_f = lc.flatten()
        mjd_time_flattened = lc_f.time.mjd - lc_f.time[0].mjd
        print(target_identifier, 'sector', sector_number, 'lightcurve flattened')
        lc_fbin = lc_f.bin(time_bin_size = BINNING)
        mjd_time_flattened_binned = lc_fbin.time.mjd - lc_fbin.time[0].mjd
        print(target_identifier, 'sector', sector_number, 'lightcurve binned and flattened')
        
        # To see this effect better we will bin up the light curves
        plt.figure(figsize=figure_size)
        plt.errorbar(mjd_time_flattened, lc_f.flux, lc_f.flux_err.data, fmt='.', ecolor='LightGrey', label='Unbinned Flattened Data')
        #print(target_identifier, 'sector', sector_number, 'lightcurve flattened data plotted')
        
        plt.errorbar(mjd_time_binned, lc_bin.flux, lc_bin.flux_err.data, fmt='.', ecolor='LightGrey', label='Binned Unflattened Data')
        #print(target_identifier, 'sector', sector_number, 'lightcurve binned data plotted')
        
        plt.errorbar(mjd_time_flattened_binned, lc_fbin.flux, lc_fbin.flux_err.data, fmt='.', ecolor='LightGrey', label='Binned Flattened Data')
        #print(target_identifier, 'sector', sector_number, 'lightcurve binned and flattened data plotted')
        
        plt.xlabel('Time [days]')
        plt.ylabel('Normalized Flux')
        plt.title('{0}, Sector {1}, Binning {2} days'.format(target_identifier,sector_number,BINNING))
        plt.legend()
        if y_limits != False:
            plt.ylim(y_limits[0],y_limits[1])
        plt.show()
        
        
        
