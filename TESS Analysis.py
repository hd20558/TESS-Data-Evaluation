# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 12:11:51 2022

@author: Louis Eddershaw
"""
import lightkurve as lk
import lightkurve.periodogram
import matplotlib.pyplot as plt
import numpy as np
import math

figure_size = (15,5)
BINNING = 0.05 #days
#binning_flag = False
detailed_analysis_threshold = 10
#y_limits = False #[0.95,1.05] #False
                #       1               2           3               4           5               6           7               8           9           10          11              12              13          14              15              16              17          18          19              20              21          22              23          24          25          26          27              28          29              30          31              32          33
observations = [[274039311,31],[126947245,1],[341849173,7],[441461124,2],[273985862,1],[266980320,1],[301407485,31],[9967126,4],[120896927,4],[231702397,1],[279614617,9],[103633672,14],[413248763,8],[441791294,14],[20291794,23],[321103619,20],[110996418,10],[149950920,21],[55525572,4],[207468071,23],[183596242,3],[350618622,2],[144539611,4],[11465798,9],[1030783,4],[7809321,4],[20448500,9],[23609565,9],[14336130,20],[441420236,1],[136916387,12],[169226822,9],[160708862,38]]
#observations = np.array(observations)
#observations = observations[observations[:,1].argsort()]
#desired_indexes = [index for index in range(len(observations))]
desired_indexes = [6,10,15,22,29]
          #TIC ID    start, end (days)
initial_masks = [[441420236,[0.7,   0.9,
                     1.2,   1.35,
                     1.5,   1.55,
                     2.55   ,2.65,
                     2.75,  2.85,
                     3.3,   3.5,
                     4.15,  4.225,
                     4.325, 4.55,
                     4.85,  4.9,
                     5.85,  5.95,
                     7.9,   8.0,
                     8.15,  8.25,
                     8.75,  8.225,
                     9.15,  9.55,
                     9.75,  9.85,
                     10.2,  10.35,
                     11.0,  11.2,
                     12.54, 13.75,
                     13.875,13.94,
                     19.4,  19.7,
                     20.05, 20.15,
                     20.2,  20.35,
                     20.55, 20.65,
                     24.05, 24.5]],
         [279614617,[4.8,   5.2,
                     17.4,  18.0,
                     21.05, 21.3]],
         [321103619,[20.4,  20.6]],
         [144539611,[6.35,6.6,
                     10.325,10.45,
                     10.875,11.0,
                     14.05,14.2,
                     14.3,14.5,
                     14.6,14.95,
                     17.8,18.0,
                     22.35,23.0]],
         [301407485,[]]]
                    #[0.5,   0.8,
                    #1.8,   2.175,
                    #3.15,  3.5,
                    #4.5,   4.8,
                    #5.8,   6.15,
                    #7.125, 7.45,
                    #8.425  ,8.8,
                    #9.8,   10.15,
                    #11.15, 11.5,
                    #15.15, 15.45,
                    #16.45, 16.8,
                    #17.8,  18.15,
                    #19.15, 19.45,
                    #20.45, 20.8,
                    #21.8,  22.15,
                    #23.15, 23.45,
                    #24.4,  24.8]
                    
secondary_masks = [[441420236,[[0.0,1.5,
                                4.5,6.3,
                                9.3,11.25,
                                14.25,16.1,         #masking the small peaks
                                19.1,20.9,
                                24.0,25.75],
                               [1.5,4.5,
                                6.3,9.3,
                                11.25,14.25,
                                16.1,19.1,          #masking the big peaks
                                20.9,24.0,  
                                25.75,28]]],
         [279614617,[]],
         [321103619,[]],
         [144539611,[[1.2,1.4,
                      3.175,3.35,
                      5.175,5.35,
                      7.15,7.35,
                      11.15,11.35,
                      15.15,15.3,                   #masking both eclipses
                      17.12,17.275,
                      19.1,19.3,
                      21.1,21.275,
                      23.075,23.25,
                      25.0,25.5],
                     [1.2,1.4,
                      5.175,5.35,
                      17.12,17.275,                         #masking secondary eclipses
                      21.1,21.275,
                      25.0,25.5],
                     [3.175,3.35,
                      7.15,7.35,
                      11.15,11.35,                          #masking primary eclipses
                      15.15,15.3,
                      19.1,19.3,
                      23.075,23.25]]],
         [301407485,[[0.5,0.75,
                      1.85,2.15,
                      3.175,3.45,
                      4.475,4.775,
                      5.8,6.1,
                      7.15,7.45],
                     [
                         ],
                     [
                         ]]]]

#analyse only the IDs at the given indexes
observations_to_analyse = [observations[desired_index] for desired_index in desired_indexes]
observations_to_analyse = np.array(observations_to_analyse)
#sort by sector number
observations_to_analyse = observations_to_analyse[observations_to_analyse[:,1].argsort()]

#loop through all desired observations
for observation_index in range(len(observations_to_analyse)):
    target_identifier = 'TIC {0}'.format(observations_to_analyse[observation_index][0])
    sector_number = observations_to_analyse[observation_index][1]
    print(observation_index + 1, 'out of', len(observations_to_analyse))

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
    # You can play around with the sigma value - what does this do to the light curve?
    # This also removes nans
    #print(lc.flux.shape, lc.flux[0:10])
    print(target_identifier, 'sector', sector_number, 'lightcurve outliers removed')
    
    #time relative to start of sector
    mjd_time = lc.time.mjd - lc.time[0].mjd
    
    ## MASK LIGHTCURVE
    mask_locations = []
    for mask_objects in initial_masks:
        if mask_objects[0] == observations_to_analyse[observation_index][0]:
            #get the list of mask timings for this observation
            mask_locations = mask_objects[1]
            
    #if there exists a mask for this observation, apply it
    if len(mask_locations) >=2:
        mask_label = ", Masked"
        masked_fluxes = list(lc.flux)
        masked_flux_errors = list(lc.flux_err)
        masked_times = list(lc.time.mjd)
        mask_index = -2
        
        #loop iterates from the last mask to the first so that this code can deal with removing entries from the flux list if desired
        for j in range(len(mjd_time) - 1,-1,-1):
            #for times within the mask, set the fluxes and errors to nothing, essentially removes it from the list while preserving the shape of the list
            if mjd_time[j] >= mask_locations[mask_index] and mjd_time[j] <= mask_locations[mask_index + 1]:
                lc.flux[j] = np.nan
                lc.flux_err[j] = np.nan
            #go to the next mask when the beginning of the period is reached and it is not the last (first) mask in the list
            elif mjd_time[j] < mask_locations[mask_index] and abs(mask_index) < len(mask_locations):
                mask_index -= 2
    else:
        mask_label = ", Unmasked"
    fig = plt.figure(figsize=figure_size)
    plt.errorbar(mjd_time, lc.flux, lc.flux_err.data, fmt='.', ecolor='LightGrey',label='Unbinned Data')
    
    
    lc_bin_flag = False
    if input('Bin Data? (y/[n]): ') != '':
        lc_bin_flag = True
        lc_bin = lc.bin(time_bin_size = BINNING) # Default time is in days
        mjd_time_binned = lc_bin.time.mjd - lc_bin.time[0].mjd
        print(target_identifier, 'sector', sector_number, 'lightcurve data binned')
        plt.errorbar(mjd_time_binned, lc_bin.flux, lc_bin.flux_err.data, fmt='.', ecolor='green',label='Binned Data')
        #print(target_identifier, 'sector', sector_number, 'lightcurve binned data plotted')
        plt.title('{0}, Sector {1}, Binning {2} days{3}'.format(target_identifier,sector_number,BINNING,mask_label))
        plt.legend()
    else:
        plt.title('{0}, Sector {1}{2}'.format(target_identifier,sector_number,mask_label))
    plt.xlabel('Time [days]')
    plt.ylabel('Normalized Flux')
    #if y_limits != False:
    #    plt.ylim(y_limits[0],y_limits[1])
    plt.show()
    
    
    
    
    
    
    
    # We can calculate the standard deviation of the data and compare it to the uncertainty on the flux given
    #scatter = np.std(lc.flux) 
    #print('Standard deviation = ', scatter)
    #print('Mean err = ', np.mean(lc.flux_err))
    
    
    
    
    #plot the errors on the flux with time rather than the actual flux value, kinda like a residual?
    if input('Plot errors on flux? (y/[n]): ') != '':
        fig = plt.figure(figsize=figure_size)
        plt.plot(mjd_time_binned,lc_bin.flux_err)
        plt.xlabel('Time [days]')
        plt.ylabel('Error on Normalised Flux Value')
        plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
        #plt.savefig('{0}-Sector-{1}'.format(target_identifier,sector_number))
        plt.show()
    
    
    
    if input('Plot Lomb-Scargle Periodogram? (y/[n]): ') != '':
        fig = plt.figure(figsize=figure_size)
        if lc_bin_flag != False:
            ls_pg = lk.periodogram.LombScarglePeriodogram.from_lightcurve(lc_bin)
            plt.plot(ls_pg.period, ls_pg.power, label='Binned Data')
        else:
            ls_pg = lk.periodogram.LombScarglePeriodogram.from_lightcurve(lc)
            plt.plot(ls_pg.period, ls_pg.power, label='Unbinned Data')
            
        ls_max_power_period = ls_pg.period_at_max_power
        print(target_identifier, 'sector', sector_number, 'Lomb-Scargle periodogram generated')
        
        plt.xlabel('Period [Days]')
        plt.xlim(0,10)
        plt.ylabel('Lomb-Scargle Power')
        plt.title('{0}, Sector {1} Lomb-Scargle Periodogram'.format(target_identifier,sector_number))
        plt.legend()
        plt.show()
    
    
    if input('Plot BLS Periodogram? (y/[n]): ') != '':
        fig = plt.figure(figsize=figure_size)
        if lc_bin_flag != False:
            bls_pg = lk.periodogram.BoxLeastSquaresPeriodogram.from_lightcurve(lc_bin)
            plt.plot(bls_pg.period, bls_pg.power, label='Binned Data')
        else:
            bls_pg = lk.periodogram.BoxLeastSquaresPeriodogram.from_lightcurve(lc)
            plt.plot(bls_pg.period, bls_pg.power, label='Unbinned Data')
        bls_max_power_period = bls_pg.period_at_max_power
        print(target_identifier, 'sector', sector_number, 'Box Least Squares periodogram generated')
        
        plt.xlabel('Period [Days]')
        plt.ylabel('BLS Power')
        plt.title('{0}, Sector {1} Box Least Squares Periodogram'.format(target_identifier,sector_number))
        plt.legend()
        plt.show()
    
    # if input('Plot Fourier Transform? (y/[n]): ') != '':
    #     sampling_rate_seconds = 1/120 #Hz, s^-1
    #     sampling_rate_days = sampling_rate_seconds * 3600 * 24 #day^-1
        
    #     flux_copy = []
    #     for i in range(len(lc.flux)):
    #         if not math.isnan(lc.flux[i]):
    #             flux_copy.append(lc.flux[i])
    #         else:
    #             flux_copy.append(0)
                
    #     f = abs(np.fft.fft(flux_copy))
    #     #f = f[:len(f)//2]
    #     freq = np.fft.fftfreq(np.array(flux_copy).shape[-1])
    #     #freq = freq[:len(freq)//2]
    #     period = 1 / freq
    #     spectrum = f.real*f.real + f.imag*f.imag
    #     nspectrum = spectrum / spectrum[0]
        
    #     plt.plot(period, nspectrum)
    #     #plt.xlim(0,30)
    #     plt.show()
    
    
        # ft = rfft(lc.flux)
        # fluxes = lc.flux.value
        # new_fluxes = []
        # for flux in fluxes:
        #     new_fluxes.append(flux)
        # #print(lc.flux)
        # duration = (lc.time[-1].mjd - lc.time[0].mjd) * 3600 * 24
        # sample_rate = len(lc.flux) / duration
        # quit()


    if input('Fold lightcurve about strongest frequency? (y/[n]): ') != '':
        print(bls_max_power_period)
        #fig = plt.figure(figsize=figure_size)
        lc.fold(bls_max_power_period).scatter(label='Period: {0:.4g}'.format(bls_max_power_period))
        plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
        plt.show()
        
        # fig = plt.figure(figsize=figure_size)
        # plt.errorbar(lc_fold.time.value, lc_fold.flux, lc_fold.flux_err.data, fmt='.', ecolor='green', label='Period: {0:.4g}'.format(bls_max_power_period))
        # plt.xlabel('Phase [JD]')
        # plt.ylabel('Normalised Flux')
        # plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
        # plt.legend()
        # plt.show() 
    
    
    if input('Inspect each day? (y/[n]): ') != '':
        days = math.floor(lc.time[-1].mjd - lc.time[0].mjd)
        day_counter = 0
        previous_day_end_index = 0
        for i in range(len(lc.time.mjd)):
            current_time = lc.time[i].mjd - lc.time[0].mjd
            if current_time >= day_counter + 1:
                plt.scatter(lc.time[previous_day_end_index:i].mjd - lc.time[0].mjd,lc.flux[previous_day_end_index:i])
                plt.title('{0}, Sector {1}, Day {2}-{3}'.format(target_identifier,sector_number,i,i+1))
                plt.xlabel('Time [days]')
                plt.ylabel('Normalized Flux')
                plt.show()
                day_counter += 1
                previous_day_end_index = i

    
    
    if input('Flatten light curve? (y/[n]): ') != '':
        # Let's define a new lightcurve ojbect where we have flattened the light curve
        lc_f = lc.flatten()
        mjd_time_flattened = lc_f.time.mjd - lc_f.time[0].mjd
        print(target_identifier, 'sector', sector_number, 'lightcurve flattened')
        
        plt.figure(figsize=figure_size)
        plt.errorbar(mjd_time_flattened, lc_f.flux, lc_f.flux_err.data, fmt='.', ecolor='LightGrey', label='Unbinned Flattened Data')
        
        if input('Bin flattened light curve? (y/[n]): ') != '':
            lc_fbin = lc_f.bin(time_bin_size = BINNING)
            mjd_time_flattened_binned = lc_fbin.time.mjd - lc_fbin.time[0].mjd
            print(target_identifier, 'sector', sector_number, 'lightcurve binned and flattened')
            
            plt.errorbar(mjd_time_flattened_binned, lc_fbin.flux, lc_fbin.flux_err.data, fmt='.', ecolor='LightGrey', label='Binned Flattened Data\n{0} d'.format(BINNING))

        
        plt.xlabel('Time [days]')
        plt.ylabel('Normalized Flux')
        plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
        plt.legend()
        #if y_limits != False:
        #    plt.ylim(y_limits[0],y_limits[1])
        plt.show()
        
        
        if input('Sequentially bin the lightcurve? (y/[n]): ') != '':
            total_data_points = len(lc_f.flux)
            total_duration = lc_f.time[-1].mjd - lc_f.time[0].mjd
            max_number_of_bins = total_data_points // 100
            #mean_errors = []
            
            binnings = []
            shot_noises = []
            bin_errors = []
            BASE = 1.5
            min_power = math.log(0.05, BASE)
            max_power = math.log(total_duration, BASE)
            for sequential_binning in np.logspace(min_power,max_power,10,base=BASE):#np.linspace(1, max_number_of_bins,10):
                lc_f_seq_bin = lc_f.bin(time_bin_size = sequential_binning)
                no_of_bins = total_duration / sequential_binning
                data_points_per_bin = total_data_points / no_of_bins
                lc_f_seq_bin.flux[-1] = np.nan
                lc_f_seq_bin.flux_err[-1] = np.nan
                #lc_f_seq_bin.flux = lc_f_seq_bin.flux[:-1]
                #lc_f_seq_bin.flux_err = lc_f_seq_bin.flux_err[:-1]
                
                #mjd_time_flattened_seq_binned = lc_f_seq_bin.time.mjd - lc_f_seq_bin.time[0].mjd
                
                #plt.errorbar(mjd_time_flattened_seq_binned, lc_f_seq_bin.flux, lc_f_seq_bin.flux_err.data, fmt='.', ecolor='LightGrey', label='Binned Flattened Data\n{0} d'.format(round(sequential_binning,2)))
                #plt.xlabel('Time [days]')
                #plt.ylabel('Normalized Flux')
                #plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
                #plt.legend()
                #plt.show()
                
                
                binnings.append(sequential_binning)
                
                mean_flux = np.nanmean(lc_f_seq_bin.flux)
                mean_error = np.nanmean(lc_f_seq_bin.flux_err)
                #std_flux = np.nanstd(lc_f_seq_bin.flux)
                #std_error = np.nanstd(lc_f_seq_bin.flux_err)

                shot_noise = 1 / math.sqrt(data_points_per_bin)
                
                bin_errors.append(mean_error)
                shot_noises.append(shot_noise)        
            
                print('Bin Size: ', sequential_binning)
                print('Average Error on Bins: ', mean_error)
                print('Expected Shot Noise: ', shot_noise,'\n ')
                
                
            fig = plt.figure(figsize=figure_size)
            plt.scatter(binnings, shot_noises)
            plt.xlabel('Binning [days]')
            plt.ylabel('Expected Shot Noise')
            plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
            #plt.xscale('log')
            #plt.yscale('log')
            plt.show()
            
            fig = plt.figure(figsize=figure_size)
            plt.scatter(binnings, bin_errors)
            plt.xlabel('Binning [days]')
            plt.ylabel('Error on Each Bin')
            plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
            #plt.xscale('log')
            #plt.yscale('log')
            plt.show()
            
            fig = plt.figure(figsize=figure_size)
            plt.scatter(binnings, bin_errors)
            plt.scatter(binnings, shot_noises)
            plt.xlabel('Binning [days]')
            plt.ylabel('Error on Each Bin')
            plt.title('{0}, Sector {1}'.format(target_identifier,sector_number))
            plt.xscale('log')
            plt.yscale('log')
            plt.show()
            
            
    
    
    
    #mask the mask
    mask_locations = []
    for mask_objects in secondary_masks:
        if mask_objects[0] == observations_to_analyse[observation_index][0]:
            #get the list of mask timings for this observation
            mask_lists = mask_objects[1]
    
    one_mask_lightcurve_flux = lc.flux
    print(one_mask_lightcurve_flux)
    one_mask_lightcurve_flux_err = lc.flux_err
    
    for mask_locations in mask_lists:
        # remove the features in the one mask lightcurve for this current list of masking locations
        #generate Lomb-Scargle periodogram for the newly masked light curve
        #fold newly masked lightcurve by Lomb-Scargle max power period
        
        if len(mask_locations) >=2:
            mask_label = ", Masked"
            masked_fluxes = list(lc.flux)
            masked_flux_errors = list(lc.flux_err)
            masked_times = list(lc.time.mjd)
            mask_index = -2
            
            #loop iterates from the last mask to the first so that this code can deal with removing entries from the flux list if desired
            for j in range(len(mjd_time) - 1,-1,-1):
                #for times within the mask, set the fluxes and errors to nothing, essentially removes it from the list while preserving the shape of the list
                if mjd_time[j] >= mask_locations[mask_index] and mjd_time[j] <= mask_locations[mask_index + 1]:
                    lc.flux[j] = np.nan
                    lc.flux_err[j] = np.nan
                #go to the next mask when the beginning of the period is reached and it is not the last (first) mask in the list
                elif mjd_time[j] < mask_locations[mask_index] and abs(mask_index) < len(mask_locations):
                    mask_index -= 2
        else:
            mask_label = ", Unmasked"
        fig = plt.figure(figsize=figure_size)
        plt.errorbar(mjd_time, lc.flux, lc.flux_err.data, fmt='.', ecolor='LightGrey',label='Unbinned Data')
        plt.show()
        
        lc.replace_column('flux', one_mask_lightcurve_flux)
        print(lc.flux)
        print(one_mask_lightcurve_flux == lc.flux)
        lc.replace_column('flux_err', one_mask_lightcurve_flux_err)
        
        
        
        
        
        
        
