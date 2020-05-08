#!/usr/bin/env python
"""
Plot remaining hotspots from the satellite search
"""
__author__ = "Sid Mau"

# Python libraries
import os
import glob
import sys
import yaml
import numpy as np
import scipy.ndimage
import matplotlib;matplotlib.use('Agg')
#from matplotlib import rc;rc('text', usetex=True);rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import fitsio as fits
import healpy as hp
import urllib
#from skimage import io

from simple.filters import quality_filter, star_filter, galaxy_filter, color_filter, dered_mag
from simple.simple_utils import cut_isochrone_path

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector
import ugali.utils.plotting
import ugali.candidate.associate
from ugali.isochrone import factory as isochrone_factory


with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    isoname = cfg[survey]['isoname']
    isosurvey = cfg[survey]['isosurvey']
    mag_max = cfg[survey]['mag_max']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

    mode = cfg[survey]['mode']
    sim_population = cfg[survey]['sim_population']
    
    band_1 = cfg[survey]['band_1']
    band_2 = cfg[survey]['band_2']
    mag = cfg[survey]['mag']
    mag_err = cfg[survey]['mag_err']
    mag_dered = cfg[survey]['mag_dered']

save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

# construct mags
mag_1 = mag.format(band_1.upper())
mag_2 = mag.format(band_2.upper())
mag_err_1 = mag_err.format(band_1.upper())
mag_err_2 = mag_err.format(band_2.upper())
mag_dered_1 = mag_dered.format(band_1.upper())
mag_dered_2 = mag_dered.format(band_2.upper())

#------------------------------------------------------------------------------

try:
    ra, dec, mod, sig, mc_source_id, field_density, = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]) 
except:
    try:
        targ_ra, targ_dec, mod, sig, mc_source_id, = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
        field_density = None
    except:
        sys.exit('ERROR! Coordinates not given in correct format.')

#------------------------------------------------------------------------------

# Set colors
cmap_gray = matplotlib.cm.Greys
cmap_gray_mask = ListedColormap(cmap_gray(np.linspace(0.1, 1.0, 100)))
cmap_gray_mask.set_bad('white')

props = dict(facecolor='white', edgecolor='black', linewidth=1) #dict(facecolor='white', edgecolor='none', alpha=0.7)

#--------------------------------------------------------------------------

# Select target region
pix_nside_select = ugali.utils.healpix.angToPix(nside, ra, dec)
pix_nside_neighbors = np.concatenate([[pix_nside_select], hp.get_all_neighbours(nside, pix_nside_select)])

# Construct data array
data_array = []
for pix_nside in pix_nside_neighbors:
    inlist = glob.glob('{}/*_{:05d}.fits'.format(datadir, pix_nside))
    for infile in inlist:
        if not os.path.exists(infile):
            continue
        data_array.append(fits.read(infile))

data = np.concatenate(data_array)

#----------------------------------------------------------------------

# quality selection
quality_cut = quality_filter(survey, data)
data = data[quality_cut]

data = dered_mag(survey, data)

color = data[mag_dered_1] - data[mag_dered_2]
mag   = data[mag_dered_1]

stars = star_filter(survey, data)
galaxies = galaxy_filter(survey, data)

iso = isochrone_factory(name='Bressan2012', survey='des', age=12.0, z=0.0001, distance_modulus=mod, band_1='g', band_2='r')
#iso_sep = iso.separation(data[mag_1], data[mag_2])
iso_filter = cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)

# projection of image
proj = ugali.utils.projector.Projector(ra, dec)
x, y = proj.sphereToImage(data['RA'], data['DEC'])

## Filters
extension = 0.05
r0 = 3.0 * extension # 3.0
r1 = 5.0 * extension # 5.0
r2 = np.sqrt(r0**2 + r1**2)
angsep = ugali.utils.projector.angsep(ra, dec, data['RA'], data['DEC'])
inner = (angsep < r0)
outer = ((angsep > r1) & (angsep < r2))
background = (angsep > r2)

#----------------------------------------------------------------------

fig, axs = plt.subplots(1, 3, figsize=(12, 3))
fig.subplots_adjust(wspace=0.5)
#fig.suptitle(r'{} (RA, Dec, $m-M$) = ({:0.2f}, {:0.2f}, {:0.0f})'.format(name, ra, dec, mod))
#fig.suptitle('PS1 {}'.format(name))

## Stellar scatter
#ax = axs[0][0]
#plt.sca(ax)
#
#ax.scatter(x[stars & iso_filter], y[stars & iso_filter], s=1, c='k', edgecolor='none', rasterized=False)
##ax.scatter(x[stars & iso_filter & outer], y[stars & iso_filter & outer], s=1, c='r', edgecolor='none', rasterized=False)
##ax.scatter(x[stars & iso_filter & inner], y[stars & iso_filter & inner], s=1, c='c', edgecolor='none', rasterized=False)
#
##theta = np.linspace(0, 2*np.pi, 1000, endpoint=True)
##
##xs_outer = np.outer([r1, r2], np.cos(theta))
##ys_outer = np.outer([r1, r2], np.sin(theta))
##xs_outer[1,:] = xs_outer[1,::-1]
##ys_outer[1,:] = ys_outer[1,::-1]
##
##xs_inner = np.outer([0, r0], np.cos(theta))
##ys_inner = np.outer([0, r0], np.sin(theta))
##xs_inner[1,:] = xs_inner[1,::-1]
##ys_inner[1,:] = ys_inner[1,::-1]
##
##ax.fill(np.ravel(xs_outer), np.ravel(ys_outer), facecolor='r', edgecolor='none', alpha=0.3)
##ax.fill(np.ravel(xs_inner), np.ravel(ys_inner), facecolor='c', edgecolor='none', alpha=0.3)
#
#ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
#
#ax.set_xlim(0.5, -0.5)
#ax.set_ylim(-0.5, 0.5)
#ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
#ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
#
## Galactic scatter
#ax = axs[0][1]
#plt.sca(ax)
#
#ax.scatter(x[galaxies & iso_filter], y[galaxies & iso_filter], s=1, c='k', edgecolor='none', rasterized=False)
#
#ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
#
#ax.set_xlim(0.5, -0.5)
#ax.set_ylim(-0.5, 0.5)
#ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
#ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
#
## CMD
#ax = axs[0][2]
#plt.sca(ax)
#
#ax.scatter(color[stars & background], mag[stars & background], s=0.1, c='k', alpha=0.1, edgecolor='none', rasterized=True)
#ax.scatter(color[stars & outer], mag[stars & outer], s=3, c='r', edgecolor='none', rasterized=False)
#ax.scatter(color[stars & inner], mag[stars & inner], s=3, c='c', edgecolor='none', rasterized=False)
#ax.scatter(color[stars & inner & iso_filter], mag[stars & inner & iso_filter], s=5, c='c', edgecolor='none', rasterized=False)
#ugali.utils.plotting.drawIsochrone(iso, lw=1, color='cyan')
#
#ax.set_xlim(-0.3, 1.0)
#ax.set_ylim(24.0, 16.0)
#ax.set_xlabel(r'$g - r$ (mag)')
#ax.set_ylabel(r'$g$ (mag)')

# Stellar histogram
ax = axs[0]
plt.sca(ax)

bound = 0.5
steps = 100.
bins = np.linspace(-bound, bound, steps)
signal = np.histogram2d(x[stars & iso_filter], y[stars& iso_filter], bins=[bins, bins])[0]
sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
pc = ax.pcolormesh(bins, bins, convolution, cmap=cmap_gray, rasterized=True)

# search kernel
#x, y = ra_proj, dec_proj
#delta_x = 0.01
#area = delta_x**2
#smoothing = 2. / 60. # Was 3 arcmin
#bins = np.arange(-8., 8. + 1.e-10, delta_x)
#centers = 0.5 * (bins[0: -1] + bins[1:])
#yy, xx = np.meshgrid(centers, centers)
#h = np.histogram2d(x, y, bins=[bins, bins])[0]
#h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
#pc = ax.pcolormesh(bins, bins, h_g.T, cmap=cmap_gray, rasterized=True)

ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)

ax.set_xlim(0.5, -0.5)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0)
cb = fig.colorbar(pc, cax=cax)
#cb.set_label('Counts')

# Galactic histogram
ax = axs[1]
plt.sca(ax)

bound = 0.5
steps = 100.
bins = np.linspace(-bound, bound, steps)
signal = np.histogram2d(x[galaxies & iso_filter], y[galaxies & iso_filter], bins=[bins, bins])[0]
sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
pc = ax.pcolormesh(bins, bins, convolution, cmap=cmap_gray, rasterized=True)

# search kernel
#x, y = ra_proj, dec_proj
#delta_x = 0.01
#area = delta_x**2
#smoothing = 2. / 60. # Was 3 arcmin
#bins = np.arange(-8., 8. + 1.e-10, delta_x)
#centers = 0.5 * (bins[0: -1] + bins[1:])
#yy, xx = np.meshgrid(centers, centers)
#h = np.histogram2d(x, y, bins=[bins, bins])[0]
#h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
#pc = ax.pcolormesh(bins, bins, h_g.T, cmap=cmap_gray, rasterized=True)

ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)

ax.set_xlim(0.5, -0.5)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0)
cb = fig.colorbar(pc, cax=cax)
#cb.set_label('Counts')

# Hess
ax = axs[2]
plt.sca(ax)

xbins = np.arange(-0.3, 1.1, 0.1)
ybins = np.arange(16., 24.0 + 0.5, 0.5)
foreground = np.histogram2d(color[stars & inner], mag[stars & inner], bins=[xbins, ybins])
background = np.histogram2d(color[stars & outer], mag[stars & outer], bins=[xbins, ybins])
fg = foreground[0].T
bg = background[0].T
fg_abs = np.absolute(fg)
bg_abs = np.absolute(bg)
mask_abs = fg_abs + bg_abs
mask_abs[mask_abs == 0.] = np.nan # mask common zeroes
signal = fg - bg
signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan
pc = ax.pcolormesh(xbins, ybins, signal, cmap=cmap_gray_mask, rasterized=True)
ugali.utils.plotting.drawIsochrone(iso, lw=1, color='cyan')

ax.set_xlim(-0.3, 1.0)
ax.set_ylim(24.0, 16.0)
ax.set_xlabel(r'$g - r$ (mag)')
ax.set_ylabel(r'$g$ (mag)')

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0)
cb = fig.colorbar(pc, cax=cax)
#cb.set_label('Counts')

#fig.savefig('./{}.pdf'.format(name), bbox_inches='tight')
file_name = 'candidate_{:0.2f}_{:0.2f}'.format(ra, dec)
fig.savefig('{}/{}.png'.format(save_dir, file_name), bbox_inches='tight')
plt.close(fig)
