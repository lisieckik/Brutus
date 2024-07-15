import h5py, gc, time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
from matplotlib import colormaps
from scipy.stats import binned_statistic_2d
from astropy.io import fits
from astropy.table import Table
import sys
import pathlib
import os


# my colormaps
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define the blue tones colormap
blue_cmap = mcolors.LinearSegmentedColormap.from_list("",[(0.0, 0.0, 0.0),
                                                          (0.05, 0.1, 0.2),
                                                          (0.2, 0.41568627450980394, 0.7764705882352941),
                                                          (0.4117647058823529, 0.5764705882352941, 0.8470588235294118),
                                                          (0.615686274509804, 0.7254901960784313, 0.9019607843137255),
                                                          (1.0, 1.0, 1.0)])

# Define the red tones colormap
red_cmap = mcolors.LinearSegmentedColormap.from_list("",
                                                     [(0.0, 0.0, 0.0),
                                                      (0.2, 0.1, 0.1),
                                                      (0.4, 0.2, 0.2),
                                                      (0.6, 0.4, 0.4),
                                                      (0.8, 0.6, 0.6),
                                                      (1.0, 0.8, 0.8),
                                                      (1.0, 1.0, 1.0)])


#cmap_rev = cmap.reversed()
#cmap2_rev = cmap2.reversed()


def rotation_matrix(direction, alpha):
    # depending on direction it gives back rotation matrix
    if direction == 'x':
        return np.matrix([[1, 0, 0],
                          [0, np.cos(alpha), -np.sin(alpha)],
                          [0, np.sin(alpha), np.cos(alpha)]])
    elif direction == 'y':
        return np.matrix([[np.cos(alpha), 0, np.sin(alpha)],
                          [0, 1, 0],
                          [-np.sin(alpha), 0, np.cos(alpha)]])
    elif direction == 'z':
        return np.matrix([[np.cos(alpha), -np.sin(alpha), 0],
                          [np.sin(alpha), np.cos(alpha), 0],
                          [0, 0, 1]])
    else:
        raise ValueError('It is not a valid direction!')


def rotation(data, alpha=0, beta=0, gamma=0):
    # uses rotation matrix to rotate all points
    try:
        values = data.values
    except:
        values = np.array(data)
    if alpha != 0:
        alpha = np.radians(alpha)
        rotm = rotation_matrix('z', alpha)
        values = np.transpose(values)
        values = np.dot(rotm, values)
        values = np.transpose(values)

    if beta != 0:
        beta = np.radians(beta)
        rotm = rotation_matrix('y', beta)
        values = np.transpose(values)
        values = np.dot(rotm, values)
        values = np.transpose(values)

    if gamma != 0:
        gamma = np.radians(gamma)
        rotm = rotation_matrix('x', gamma)
        values = np.transpose(values)
        values = np.dot(rotm, values)
        values = np.transpose(values)

    return pd.DataFrame(values, columns=['x', 'y', 'z'])


def segmentation(catalog, alpha, beta, gamma, segment_size=1e7, canvas_size=1500, p_type='particle1',
                 weight=None, x0=6000, xk=16000, y0=6000, yk=16000, z0=6000, zk=16000):
    # divides dataset in segments, limits it around your galaxy, rotates it and make 2d histogram

    # changes the frame of reference, that the galaxy is in the center
    # it makes rotation easier
    xmid = (xk - x0) / 2
    ymid = (yk - y0) / 2
    zmid = (zk - z0) / 2

    # finding all vertexes to rotate them and find x and y ranges for plot
    cube_range = np.array([[-xmid, -ymid, -zmid], [xmid, -ymid, -zmid],
                           [-xmid, ymid, -zmid], [xmid, ymid, -zmid],
                           [-xmid, -ymid, zmid], [xmid, -ymid, zmid],
                           [-xmid, ymid, zmid], [xmid, ymid, zmid]])
    cube_range = pd.DataFrame(cube_range, columns=['x', 'y', 'z'])
    cube_range = rotation(cube_range, alpha=alpha, beta=beta, gamma=gamma)

    x = np.linspace(np.min(cube_range['x']), np.max(cube_range['x']), canvas_size + 1)
    y = np.linspace(np.min(cube_range['y']), np.max(cube_range['y']), canvas_size + 1)

    # dividing dataset in segments
    segment_size = int(segment_size)
    n_iter = int(len(catalog[p_type]) / segment_size) + 1
    print('Number of particles=', len(catalog[p_type]), n_iter)
    # building empty array for resutls
    z = np.zeros([canvas_size, canvas_size])
    for j in range(n_iter):
        # finding ranges for this segment
        amin = int(j * segment_size)
        if j + 1 == n_iter:
            amax = len(catalog[p_type])
            print(j+1)
        else:
            amax = int((j + 1) * segment_size)

        # taking only particles in the range
        data_bari = pd.DataFrame(catalog[p_type][amin:amax], columns=['x', 'y', 'z'])

        # taking only particles in this galaxy
        data_bari = data_bari[(data_bari['x']>=x0) & (data_bari['x']<=xk)]
        data_bari = data_bari[(data_bari['y']>=y0) & (data_bari['y']<=yk)]
        data_bari = data_bari[(data_bari['z']>=z0-500) & (data_bari['z']<=zk+500)]

        # chanching their frame of reference
        data_bari['x'] = data_bari['x'] - (xmid + x0)
        data_bari['y'] = data_bari['y'] - (ymid + y0)
        data_bari['z'] = data_bari['z'] - (zmid + z0)

        # rotating
        rot_data_bari = rotation(data_bari, alpha=alpha, beta=beta, gamma=gamma)


        # looking for weights
        if weight == None:
            w = [1] * len(rot_data_bari)
        else:
            print('amax', amax)
            w = weight[amin:amax][data_bari.index]
            print(len(w))
            print(len(rot_data_bari))
            # making a 2d hist
        grid, xx, yy, _ = binned_statistic_2d(rot_data_bari['x'].values,
                                              rot_data_bari['y'].values, w, 'sum',
                                              bins=[canvas_size, canvas_size], range=[[x[0], x[-1]], [y[0], y[-1]]])

        z += grid

    return y, x, z


# this depends on the RAM memory, 1e8 is enough
segment_size = 1e8

# rotation angles
beta = 0 # around y direction (long screen side)
alpha = 60 #around z_axis (short screen side)
gamma = 30 #around x direction (through the screen)
# 10 is optimal
cs = 10
# dx is physical size of the cube in similuation
# fsize is size of the figure, the best option is to keep canvas_size 10*fsize
# position of your galaxy from catalog (*68 to go from comoving to physical)
pos = np.array([17957.41,  22369.352, 27080.266]) * 0.68



for dx, fsize in zip([80], [30]):
    x, y, z = pos
    # how many bins will be there (like 2dhist)
    canvas_size = fsize * cs

    # centering on the galaxy
    fig_ranges = [[x - dx, x + dx], [y - dx, y + dx], [z - dx, z + dx]]

    # snapshot ID
    n = 129
    # file from simba
    redshifts = np.loadtxt('plik.txt')[:, 1]
    # catalog from simba
    file = "/home/lorenzong/Desktop/imaging_simba/high_res_box/snap_m25n512_%i.hdf5" % n
    #file = "/media/lorenzong/Data1/simba_hig_res/snap_m25n512_%i.hdf5" % n
    catalog = h5py.File(file, 'r')

    # part 0 -- gas
    # part 1 -- dark matter
    # part 4 -- stars
    # part 5 -- black hole
    cmaps = [colormaps["inferno"], colormaps["viridis"]]
    ptypes = ['PartType0', 'PartType4']
    titles = ['Gas component', 'Stellar component']
    fig, ax = plt.subplots(1, 2, figsize=(3.5*fsize/10, 1.5*fsize/10), sharex=True, sharey=True)
    for p in range(len(ptypes)):
        coords = catalog[ptypes[p]]
        # magic
        x, y, z = segmentation(coords, alpha, beta, gamma,
                               canvas_size=canvas_size, segment_size=segment_size, p_type='Coordinates',
                               x0=fig_ranges[0][0], xk=fig_ranges[0][1],
                               y0=fig_ranges[1][0], yk=fig_ranges[1][1],
                               z0=fig_ranges[2][0], zk=fig_ranges[2][1])#, weight=catalog[ptypes[p]]['Sigma'])


        try:
            z_new += z #(z-z.min())/(z.max()-z.min())
        except NameError:
            z_new = z #(z-z.min())/(z.max()-z.min())
        print(z_new[np.where(z_new==np.nan)[0]])
        #z_new[z==0] = np.nan
        cmaps[0].set_bad(alpha=0)
        #ax[p].set_position([0.1, 0.1, 0.8, 0.8])
        pcm = ax[p].pcolormesh(x, y, z_new + 1, cmap=cmaps[p], norm='log')
        # Adjust tick parameters
        #ax[p].tick_params(axis='both', which='major', labelsize=80, width=10, length=50)  # Set the size of major ticks
        #ax[p].tick_params(axis='both', which='minor', labelsize=80, width=10, length=50)  # Set the size of minor ticks
        ax[p].text(0.05, 0.95, titles[p], transform=ax[p].transAxes, va='top', ha='left', color='w')
        ax[p].set_xlabel('Kpc')
        fig.tight_layout()
        if p==1:
            cbar = fig.colorbar(pcm, ax=ax[p], orientation='vertical', label='# paricles')
        else:
            cbar = fig.colorbar(pcm, ax=ax[p], orientation='vertical')
            ax[0].set_ylabel('Kpc', x=0.06)
        del z_new
    plt.savefig(f'/home/lorenzong/Desktop/imaging_simba/out/single_gal/massive_agns/ztest_image.png', transparent=False)
    plt.close()
