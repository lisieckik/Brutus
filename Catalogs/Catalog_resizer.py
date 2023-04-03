import h5py
import yt
import caesar
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import time
import psutil
import gc
import os

z_snap = np.flip(np.genfromtxt('/home/lorenzong/Desktop/Catalogs/SIMBA_catalogs/z_simba_snap', usecols=1, unpack=True))

# names of columns in result fits file (first galaxies then halo)
names = ['parent_ID', 'result_ID', 'redshift',
         'pos_x', 'pos_y','pos_z',
         'vel_x','vel_y','vel_z',
         'L_FIR', 'SFR', 'stellar_mass',
         'gas_mass','HI_mass','H2_mass',
         'dust_mass','bh_mass',
         'radii_gas_half_mass','radii_stellar_half_mass',
         'radii_total_half_mass', 'metallicity_mass_weighted',
         'metallicity_sfr_weighted', 'metallicity_stellar',
         'velocity_dispersions_gas', 'velocity_dispersions_stellar',
         'ages_mass_weighted', 'ages_metal_weighted',  'central','halo_ID',
         'halo_position_x', 'halo_position_y', 'halo_position_z',
         'halo_velocity_x', 'halo_velocity_y', 'halo_velocity_z',
         'halo_circular_velocity', 'halo_r200', 'halo_temperature',
         'halo_stellar_metallicities', 'halo_masses_HI', 'halo_masses_H2',
         'halo_masses_gas', 'halo_masses_dm', 'halo_masses_dust',
         'halo_total_mass', 'halo_bh_mass'
         ]

# names of attributes collected from the snapshot catalog
attributes_galaxies = ['pos',
                       'vel',
                       'L_FIR',
                       'sfr',
                       ['dicts', 'masses.stellar'],
                       ['dicts', 'masses.gas'],
                       ['dicts', 'masses.HI'],
                       ['dicts', 'masses.H2'],
                       ['dicts', 'masses.dust'],
                       ['dicts', 'masses.bh'],
                       ['dicts', 'radii.gas_half_mass'],
                       ['dicts', 'radii.stellar_half_mass'],
                       ['dicts', 'radii.total_half_mass'],
                       ['dicts', 'metallicities.mass_weighted'],
                       ['dicts', 'metallicities.sfr_weighted'],
                       ['dicts', 'metallicities.stellar'],
                       ['dicts', 'velocity_dispersions.gas'],
                       ['dicts', 'velocity_dispersions.stellar'],
                       ['dicts','ages.mass_weighted'],
                       ['dicts','ages.metal_weighted'],
                       'central',
                       ]


proprety_list_halo = ['GroupID',
                 'minpotpos',
                 'minpotvel',
                 ['dicts', 'virial_quantities.circular_velocity'],
                 ['dicts', 'virial_quantities.r200'],
                 ['dicts', 'virial_quantities.temperature'],
                 ['dicts', 'metallicities.stellar'],
                 ['dicts', 'masses.HI'],
                 ['dicts', 'masses.H2'],
                 ['dicts', 'masses.gas'],
                 ['dicts', 'masses.dm'],
                 ['dicts', 'masses.dust'],
                 ['dicts', 'masses.total'],
                 ['dicts', 'masses.bh'],
                 ]

# calculating the number of all needed columns in result array
number_of_columns = 3
for i in range(len(attributes_galaxies)):
    if attributes_galaxies[i] == 'pos' or attributes_galaxies[i] == 'vel':
        number_of_columns+=3
    else:
        number_of_columns+=1

for i in range(len(proprety_list_halo)):
    if proprety_list_halo[i] == 'minpotpos' or proprety_list_halo[i] == 'minpotvel':
        number_of_columns+=3
    else:
        number_of_columns+=1
        
#read file about protoclusters
protocluster_cat = '/home/lorenzong/Desktop/Catalogs/protoclusters_z2.hdf5'
#protocluster_cat = './field_sample.hdf5'
protocluster_cat = h5py.File(protocluster_cat, 'r')

# making empty header HDU for fits
hdr = fits.Header()
empty_primary = fits.PrimaryHDU(header=hdr)

for i in range(75):
    ind_snap = str(151 - i) # write the name of the snapshot

    if len(ind_snap) < 3:
        ind_snap = '0' * (3 - len(ind_snap)) + ind_snap # if snapshot number is less than 100 we need to add 0 in front
    print(ind_snap)

    # IDs of galaxies at z=0
    protocluster_IDs = list(protocluster_cat[ind_snap].keys())

    # open corresponding simba catalog catalog 
    file = '/home/lorenzong/Desktop/Catalogs/SIMBA_catalogs/m100n1024_%s.hdf5' % ind_snap
    sim = h5py.File(file)

    # getting the redshift information
    redshift = round(z_snap[i], 3)

    # collecting information about number of ALL the progenitors in the snapshot
    number_of_IDs = len(protocluster_IDs)
    all_ID_array = []
    progenitors_ID = []
    for j in range(number_of_IDs):
        #print(protocluster_IDs[j], protocluster_cat[ind_snap][protocluster_IDs[j]])
        IDs = protocluster_cat[ind_snap][protocluster_IDs[j]][:]
        all_ID_array += list(IDs)
        progenitors_ID += len(IDs[IDs>0]) * [protocluster_IDs[j]]
    all_ID_array = np.array(all_ID_array).astype(int)
    all_ID_array = all_ID_array[all_ID_array>0]
    progenitors_ID = np.array(progenitors_ID).astype(int)
    #progenitors_ID = progenitors_ID[progenitors_ID>0]
    row_size = len(all_ID_array)
    

    # new array for all the results
    new_cat = np.empty([row_size,number_of_columns])
    row_size = 0

    # collecting all the data from SIMBA snapshot catalog
    new_cat[:, 0] = all_ID_array
    new_cat[:, 1] = progenitors_ID
    new_cat[:, 2] = len(all_ID_array)*[redshift]
    print('=================', np.shape(new_cat))
    
    additional_column = 3
    matchgal = all_ID_array[all_ID_array>0]
    for k in range(len(attributes_galaxies)):
        if len(attributes_galaxies[k]) ==2:
            new_cat[:, k + additional_column] =\
                np.array(sim['galaxy_data'][attributes_galaxies[k][0]][attributes_galaxies[k][1]])[matchgal]
        else:
            if attributes_galaxies[k] == 'pos' or attributes_galaxies[k] == 'vel':
                print('matchgal', np.shape(matchgal))
                print('right', np.shape(np.array(sim['galaxy_data'][attributes_galaxies[k]])[matchgal]))
                print('left', np.shape(new_cat[:, k + additional_column: k + additional_column +3]))
                new_cat[:, k + additional_column: k + additional_column +3] = \
                    np.array(sim['galaxy_data'][attributes_galaxies[k]])[matchgal]
                additional_column +=2
            else:
                new_cat[:, k+ additional_column] =\
                    np.array(sim['galaxy_data'][attributes_galaxies[k]])[matchgal]



    halo_IDs = np.array(sim['galaxy_data']['parent_halo_index'])[matchgal]
    additional_column += len(attributes_galaxies)
    for k in range(len(proprety_list_halo)):
        if len(proprety_list_halo[k]) == 2:
            new_cat[:, k+additional_column] = np.array(sim['halo_data'][proprety_list_halo[k][0]][proprety_list_halo[k][1]])[halo_IDs]
        else:
            if proprety_list_halo[k] == 'minpotpos' or proprety_list_halo[k] == 'minpotvel':
                new_cat[:, k + additional_column: k + additional_column + 3] = np.array(sim['halo_data'][proprety_list_halo[k]])[halo_IDs]
                additional_column += 2
            else:
                new_cat[:, k+ additional_column] = np.array(sim['halo_data'][proprety_list_halo[k]])[halo_IDs]


                
    cols = []
    for j in range(new_cat.shape[1]):
        col = fits.Column(name=names[j], array=new_cat[:,j], format='E')
        cols.append(col)

    table_hdu = fits.BinTableHDU.from_columns(cols)
    hdul = fits.HDUList([empty_primary, table_hdu])
    hdul.writeto('/home/lorenzong/Desktop/Catalogs/new_cat/%s.fits'%ind_snap, overwrite=True)
    print('RAM memory used:', psutil.virtual_memory()[2])