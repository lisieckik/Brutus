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
         'ages_mass_weighted', 'ages_metal_weighted',  'central',
         'halo_position_x', 'halo_position_y', 'halo_position_z',
         'halo_velocity_x', 'halo_velocity_y', 'halo_velocity_z',
         'halo_circular_velocity', 'halo_r200', 'halo_temperature',
         'halo_total_mass', 'halo_bh_mass'
         ]

# names of attributes collected from the snapshot catalog
attributes_galaxies = ['pos',
                       'vel',
                       'L_FIR',
                       'sfr',
                       ['masses', 'stellar'],
                       ['masses', 'gas'],
                       ['masses', 'HI'],
                       ['masses', 'H2'],
                       ['masses', 'dust'],
                       ['masses', 'bh'],
                       ['radii', 'gas_half_mass'],
                       ['radii', 'stellar_half_mass'],
                       ['radii', 'total_half_mass'],
                       ['metallicities', 'mass_weighted'],
                       ['metallicities', 'sfr_weighted'],
                       ['metallicities', 'stellar'],
                       ['velocity_dispersions', 'gas'],
                       ['velocity_dispersions', 'stellar'],
                       ['ages','mass_weighted'],
                       ['ages','metal_weighted'],
                       'central',
                       ]


proprety_list_halo = [
                 'minpotpos',
                 'minpotvel',
                 ['virial_quantities', 'circular_velocity'],
                 ['virial_quantities', 'r200'],
                 ['virial_quantities', 'temperature'],
                 ['masses', 'total'],
                 ['masses', 'bh'],
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
protocluster_cat = 'protoclusters_z2.hdf5'
protocluster_cat = h5py.File(protocluster_cat, 'r')

# making empty header HDU for fits
hdr = fits.Header()
empty_primary = fits.PrimaryHDU(header=hdr)

for i in range(4):
    pliki = os.listdir('test_new')
    nnn = 151
    for j in range(len(pliki)):
        nr = int(pliki[j].split('.')[0])
        if nr <= nnn: nnn = nr-1
    t0 = time.time()
    # ind_snap lower for snapshot catalog for collecting the data about all mergers participants
    ind_snap = str(nnn)
    if len(ind_snap) < 3:
        ind_snap = '0' * (3 - len(ind_snap)) + ind_snap

    # list of the result galaxies in higher snapshot
    protocluster_IDs = list(protocluster_cat[ind_snap].keys())

    # opening catalog of lower snapshot
    file = 'catalogs/m100n1024_%s.hdf5' % ind_snap
    sim = caesar.load(file)

    # print([getattr(j, 'ages')['mass_weighted'] for j in sim.galaxies])



    # getting the redshift information
    redshift = round(sim.simulation.redshift, 3)

    # collecting information about number of ALL the progenitors in lower snapshot
    number_of_IDs = len(protocluster_IDs)
    all_ID_array = []
    progenitors_ID = []
    for j in range(number_of_IDs):
        IDs = protocluster_cat[ind_snap][protocluster_IDs[j]][:]
        all_ID_array += list(IDs)
        progenitors_ID += len(IDs) * [protocluster_IDs[j]]
    row_size = len(all_ID_array)

    # new array for all the results
    new_cat = np.empty([row_size,number_of_columns])
    row_size = 0

    # collecting all the data from lower snapshot catalog
    new_cat[:, 0] = all_ID_array
    new_cat[:, 1] = progenitors_ID
    new_cat[:, 2] = len(all_ID_array)*[redshift]

    additional_column = 3

    for k in range(len(attributes_galaxies)):
        if len(attributes_galaxies[k]) ==2:
            new_cat[:, k + additional_column] =\
                np.array([getattr(j, attributes_galaxies[k][0])[attributes_galaxies[k][1]] for j in sim.galaxies])[all_ID_array]
        else:
            if attributes_galaxies[k] == 'pos' or attributes_galaxies[k] == 'vel':
                new_cat[:, k + additional_column: k + additional_column +3] = \
                    np.array([getattr(j, attributes_galaxies[k]) for j in sim.galaxies])[all_ID_array]
                additional_column +=2
            else:
                new_cat[:, k+ additional_column] =\
                    np.array([getattr(j,attributes_galaxies[k]) for j in sim.galaxies])[all_ID_array]



    halo_IDs = np.array([i.parent_halo_index for i in sim.galaxies])[all_ID_array]
    additional_column += len(attributes_galaxies)
    for k in range(len(proprety_list_halo)):
        if len(proprety_list_halo[k]) == 2:
            new_cat[:, k+additional_column] = np.array([getattr(j, proprety_list_halo[k][0])[proprety_list_halo[k][1]] for j in sim.halos])[halo_IDs]
        else:
            if proprety_list_halo[k] == 'minpotpos' or proprety_list_halo[k] == 'minpotvel':
                new_cat[:, k + additional_column: k + additional_column + 3] = np.array([getattr(j, proprety_list_halo[k]) for j in sim.halos])[halo_IDs]
                additional_column += 2
            else:
                new_cat[:, k+ additional_column] = np.array([getattr(j, proprety_list_halo[k]) for j in sim.halos])[halo_IDs]


    cols = []


    for j in range(new_cat.shape[1]):
        col = fits.Column(name=names[j], array=new_cat[:,j], format='E')
        cols.append(col)

    table_hdu = fits.BinTableHDU.from_columns(cols)
    hdul = fits.HDUList([empty_primary, table_hdu])
    hdul.writeto('test_new/%s.fits'%ind_snap, overwrite=True)
    print('RAM memory % used:', psutil.virtual_memory()[2])
    print('cat nr%s written'%(ind_snap), 'it took %.2f seconds' %((time.time()-t0)))

