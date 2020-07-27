#!/usr/bin/env python
"""
Diagnostic plot functions
"""
__author__ = "Sidney Mau"

import os
import yaml
import argparse
import fitsio

from astropy.coordinates import SkyCoord

import numpy as np
from scipy import interpolate
from scipy.signal import argrelextrema
import scipy.ndimage

import matplotlib.pylab as plt
import matplotlib
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
#matplotlib.use('Agg')

import simple.survey
import simple.search
import simple.isochrone

import ugali.utils.projector
import ugali.candidate.associate

################################################################################

# Set colors
cmap_gray = matplotlib.cm.Greys
cmap_color = matplotlib.cm.viridis
c_color = cmap_color(0.25)
p_color = cmap_color(1.0)
cmap_gray_mask = ListedColormap(cmap_gray(np.linspace(0.1, 1.0, 100)))
cmap_gray_mask.set_bad('white')

################################################################################

def get_iso_filter(region, data, iso):
    iso_filter = iso.cut_separation(region.survey.band_1.lower(), region.survey.band_2.lower(), data[region.survey.mag_1], data[region.survey.mag_2], data[region.survey.mag_err_1], data[region.survey.mag_err_2], radius=0.1)
    if region.survey.band_3 is not None:
        iso_filter &= iso.cut_separation(region.survey.band_2.lower(), region.survey.band_3.lower(), data[region.survey.mag_2], data[region.survey.mag_3], data[region.survey.mag_err_2], data[region.survey.mag_err_3], radius=0.1)
    return iso_filter

def get_g_radius(region, data, iso):
    """Analyze a candidate"""

    iso_filter = get_iso_filter(region, data, iso)

    # g_radius estimate
    angsep = ugali.utils.projector.angsep(region.ra, region.dec, data[region.survey.catalog['basis_1']], data[region.survey.catalog['basis_2']])

    bins = np.linspace(0, 0.4, 21) # deg
    centers = 0.5*(bins[1:] + bins[0:-1])
    area = np.pi*(bins[1:]**2 - bins[0:-1]**2) * 60**2
    hist = np.histogram(angsep[(angsep < 0.4) & iso_filter], bins=bins)[0] # counts

    f_interp = interpolate.interp1d(np.linspace(centers[0], centers[-1], len(hist)), hist/area, 'cubic')
    f_range = np.linspace(centers[0], centers[-1], 1000)
    f_val = f_interp(f_range)

    pairs = list(zip(f_range, f_val))

    peak = max(pairs[:int(len(pairs)/4)], key=lambda x: x[1]) # find peak within first quarter

    def peak_index(pairs, peak):
        for i in range(len(pairs)):
            if pairs[i] == peak:
                return i

    osc = int(0.04/0.4*1000) # +/- 0.04 (rounded down) deg oscillation about local extremum
    relmin = argrelextrema(f_val, np.less, order=osc)[0]

    try:
        if len(relmin) > 0:
            #half_point = f_range[relmin[0]]
            i = 0
            while ((f_range[relmin[i]] <= f_range[peak_index(pairs,peak)]) & (i <= len(relmin)-1)):
                i+=1
            half_point = f_range[relmin[i]]
        elif len(relmin) == 0:
            half_peak = (peak[1] + np.mean(f_val[len(f_val)/4:]))/2. # normalized to background (after first quarter)
            #half_peak = np.mean(f_val[len(f_val)/4:])
            half_pairs = []
            for i in pairs[peak_index(pairs, peak):len(pairs)/2]: # start after peak, stay within first quarter
                if i != peak:
                    half_pairs.append((i[0], abs(i[1]-half_peak)))
            half_point = min(half_pairs, key=lambda x: x[1])[0] # deg
    except:
        half_point = 0.1 # fixed value to catch errors

    g_min = 0.5/60. # deg
    g_max = 12./60. # deg

    if half_point < g_min:
        g_radius = g_min
    elif half_point > g_max:
        g_radius = g_max
    else:
        g_radius = half_point # deg

    return(g_radius)

def density_plot(ax, region, data, g_radius, iso, type):
    """Stellar density plot"""

    if type == 'stars':
        ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')
        filter = get_iso_filter(region, data, iso)
    elif type == 'galaxies':
        #ax.set_title('Galactic Density')
        ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top')
        filter = get_iso_filter(region, data, iso)
    elif type == 'blue_stars':
        filter = get_iso_filter(region, data, iso) & (data[region.survey.mag_1] - data[region.survey.mag_2] < 0.4)
        ax.text(0.05, 0.95, 'Blue stars', transform=ax.transAxes, verticalalignment='top')
    
    # projection of image
    fd = data[filter]
    ra = fd[region.survey.catalog['basis_1']]
    dec = fd[region.survey.catalog['basis_2']]
    x, y = region.proj.sphereToImage(ra, dec)

    bound = 0.5 #1.
    steps = 100.
    bins = np.linspace(-bound, bound, steps)
    signal = np.histogram2d(x, y, bins=[bins, bins])[0]
    sigma = 0.01 * (0.25 * np.arctan(0.25*g_radius*60. - 1.5) + 1.3)
    convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps))
    pc = ax.pcolormesh(bins, bins, convolution.T, cmap=cmap_gray)

    ax.set_xlim(bound, -bound)
    ax.set_ylim(-bound, bound)
    ax.set_xlabel(r'$\Delta$ RA (deg)')
    ax.set_ylabel(r'$\Delta$ Dec (deg)')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '5%', pad=0)
    plt.colorbar(pc, cax=cax, label='counts')


def star_plot(ax, region, data, iso):
    """Star bin plot"""

    iso_filter = get_iso_filter(region, data, iso)
    # projection of image
    x, y = region.proj.sphereToImage(data[iso_filter][region.survey.catalog['basis_1']], data[iso_filter][region.survey.catalog['basis_2']])

    ax.scatter(x, y, edgecolor='none', s=3, c='black')

    ax.set_xlim(0.5, -0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(r'$\Delta$ RA (deg)')
    ax.set_ylabel(r'$\Delta$ Dec (deg)')
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')

def star_plot_aperture(ax, region, data, iso, r):
    """Zoomed star bin plot, with aperture"""

    iso_filter = get_iso_filter(region, data, iso)
    # projection of image
    x, y = region.proj.sphereToImage(data[iso_filter][region.survey.catalog['basis_1']], data[iso_filter][region.survey.catalog['basis_2']])

    dot_size = np.clip(3*(.05/r), 1, 6)
    ax.scatter(x, y, edgecolor='none', s=dot_size, c='black')

    aperture = Circle(xy=(0,0), radius=r, edgecolor='blue', linewidth=1.0, linestyle = '--', fill=False, zorder=10)
    ax.add_patch(aperture)

    inner = Circle(xy=(0,0), radius=0.3, edgecolor='red', linewidth=1.0, linestyle = '--', fill=False, zorder=10)
    ax.add_patch(inner)
    outer = Circle(xy=(0,0), radius=0.5, edgecolor='red', linewidth=1.0, linestyle = '--', fill=False, zorder=10)
    ax.add_patch(outer)

    ax.set_xlim(r*2, -r*2)
    ax.set_ylim(-r*2, r*2)
    ax.set_xlabel(r'$\Delta$ RA (deg)')
    ax.set_ylabel(r'$\Delta$ Dec (deg)')
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')


def cm_plot(ax, region, data, iso, g_radius, type):
    """Color-magnitude plot"""

    if type == 'stars':
        ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')
    elif type == 'galaxies':
        ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top')

    iso_filter = get_iso_filter(region, data, iso)
    angsep = ugali.utils.projector.angsep(region.ra, region.dec, data[region.survey.catalog['basis_1']], data[region.survey.catalog['basis_2']])
    annulus = (angsep > g_radius) & (angsep < 1.)
    nbhd = (angsep < g_radius)

    # Plot background objects
    ax.scatter(data[region.survey.mag_1][annulus] - data[region.survey.mag_2][annulus], data[region.survey.mag_1][annulus], c='k', alpha=0.1, edgecolor='none', s=1)

    # Plot isochrone
    ax.plot(iso.color, iso.mag_1 + iso.distance_modulus, c='k', lw=1)

    # Plot objects in nbhd
    ax.scatter(data[region.survey.mag_1][nbhd] - data[region.survey.mag_2][nbhd], data[region.survey.mag_1][nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(g_radius))

    # Plot objects in nbhd and near isochrone
    ax.scatter(data[region.survey.mag_1][nbhd & iso_filter] - data[region.survey.mag_2][nbhd & iso_filter], data[region.survey.mag_1][nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    ax.legend(loc='upper right')

    ax.set_xlim(-0.5, 1)
    #ax.set_ylim(mag_max, 16)
    ax.set_ylim(26.0, 16)
    ax.set_xlabel('{} - {} (mag)'.format(region.survey.band_1.lower(), region.survey.band_2.lower()))
    ax.set_ylabel('{} (mag)'.format(region.survey.band_1.lower()))

def hess_plot(ax, region, data, iso, g_radius):
    """Hess plot"""

    c1 = SkyCoord(region.ra, region.dec, unit='deg')

    r_near = 2.*g_radius # annulus begins at 2*g_radius away from centroid
    r_far = np.sqrt(5.)*g_radius # annulus has same area as inner area

    angsep = ugali.utils.projector.angsep(region.ra, region.dec, data[region.survey.catalog['basis_1']], data[region.survey.catalog['basis_2']])
    inner = (angsep < g_radius)
    outer = ((angsep > r_near) & (angsep < r_far))

    xbins = np.arange(-0.5, 1.1, 0.1)
    #ybins = np.arange(16., mag_max + 0.5, 0.5)
    ybins = np.arange(16., 26.0, 0.5)
    foreground = np.histogram2d(data[region.survey.mag_1][inner] - data[region.survey.mag_2][inner], data[region.survey.mag_1][inner], bins=[xbins, ybins])
    background = np.histogram2d(data[region.survey.mag_1][outer] - data[region.survey.mag_2][outer], data[region.survey.mag_1][outer], bins=[xbins, ybins])
    fg = foreground[0].T
    bg = background[0].T
    fg_abs = np.absolute(fg)
    bg_abs = np.absolute(bg)
    mask_abs = fg_abs + bg_abs
    mask_abs[mask_abs == 0.] = np.nan # mask common zeroes
    signal = fg - bg
    signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan

    ax.set_xlim(-0.5, 1.0)
    #ax.set_ylim(mag_max, 16)
    ax.set_ylim(26.0, 16)
    ax.set_xlabel('{} - {} (mag)'.format(region.survey.band_1.lower(), region.survey.band_2.lower()))
    ax.set_ylabel('{} (mag)'.format(region.survey.band_1.lower()))

    pc = ax.pcolormesh(xbins, ybins, signal, cmap=cmap_gray_mask)
    ax.plot(iso.color, iso.mag_1 + iso.distance_modulus, lw=2, c='k', zorder=10, label='Isochrone')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '5%', pad=0)
    plt.colorbar(pc, cax=cax, label='counts')

def radial_plot(ax, region, stars, galaxies, iso, g_radius, field_density=None):
    """Radial distribution plot"""

    bins = np.linspace(0, 0.4, 21) # deg
    centers = 0.5*(bins[1:] + bins[0:-1])
    area = np.pi*(bins[1:]**2 - bins[0:-1]**2) * 60**2

    def interp_values(type, seln):
        if type == 'stars':
            data = stars
        elif type == 'galaxies':
            data = galaxies
        #data = region.get_data(type)
        iso_filter = get_iso_filter(region, data, iso)
        angsep = ugali.utils.projector.angsep(region.ra, region.dec, data[region.survey.catalog['basis_1']], data[region.survey.catalog['basis_2']])
        if seln == 'f':
            filter = iso_filter
        elif seln == 'u':
            filter = ~iso_filter

        hist = np.histogram(angsep[(angsep < 0.4) & filter], bins=bins)[0] # counts

        f_interp = interpolate.interp1d(np.linspace(centers[0], centers[-1], len(hist)), hist/area, 'cubic')
        f_range = np.linspace(centers[0], centers[-1], 1000)
        f_val = f_interp(f_range)

        val = hist/area
        yerr = np.sqrt(hist)/area

        return(f_range, f_val, val, yerr)

    f_range, f_val, val, y_err = interp_values('stars', 'f')
    pairs = list(zip(f_range, f_val))
    peak = max(pairs[:int(len(pairs)/4)], key=lambda x: x[1]) # find peak within first quarter
    def peak_index(pairs, peak):
        for i in range(len(pairs)):
            if pairs[i] == peak:
                return i
    ax.axvline(x=f_range[peak_index(pairs,peak)], color='m', label='peak')
    ax.axvline(x=g_radius, color='r', label='g_radius')

    f_range, f_val, val, yerr = interp_values('galaxies', 'f')
    ax.plot(f_range, f_val, '-g', label='Filtered Galaxies')

    f_range, f_val, val, yerr = interp_values('stars', 'u')
    ax.plot(f_range, f_val, '-k', alpha=0.25, label='Unfiltered Stars')

    f_range, f_val, val, yerr = interp_values('stars', 'f')
    ax.plot(centers, val, '.b')
    ax.errorbar(centers, val, yerr=y_err, fmt='none', ecolor='b', elinewidth=1, capsize=5)
    ax.plot(f_range, f_val, '-b', label='Filtered Stars')

    if field_density:
        ax.axhline(y=field_density, color='blue', ls='--', label='Model Field')

    ymax = plt.ylim()[1]
    ax.annotate(r'$\approx %0.1f$' + str(round(g_radius, 3)) + '$^\circ$', (g_radius*1.1, ymax/50.), color='red', bbox=dict(boxstyle='round,pad=0.0', fc='white', alpha=0.75, ec='white', lw=0))

    ax.legend(loc='upper right')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Angular Separation (deg)')
    ax.set_ylabel('Denisty (arcmin$^{-2})$')


def make_plot(survey, candidate=None, **kwargs):
    """ Creates and saves a 9-panel plot for a single candidate.
    If a candidate is passed, its parameters (coordinates, modulus, etc.)
    will be used. kwargs can be used to supercede these parameters if desired.
    If a candidate is not passed, all necessary parameters must be specified
    in kwargs.
    Parameters: sig, ra, dec, mod, r, n_obs, n_model
    """
    keys = ['sig', 'ra', 'dec', 'mod', 'r', 'n_obs', 'n_model']
    params = dict.fromkeys(keys)
    if candidate is not None:
        try: # simple
            params['sig'] = round(candidate['SIG'], 2)
        except: # ugali
            params['sig'] = round(candidate['TS'], 2)
        params['ra']      = round(candidate[survey.catalog['basis_1']], 2)
        params['dec']     = round(candidate[survey.catalog['basis_2']], 2)
        params['mod']     = round(candidate['MODULUS'], 2)
        params['r']       = round(candidate['R'], 3)
        params['n_obs']   = candidate['N_OBS']
        params['n_model'] = candidate['N_MODEL']
    for key in keys:
        try:
            params[key] = kwargs[key]
        except KeyError: # key not in kwargs
            if params[key] is None: # Missing argument
                raise TypeError('{} parameter required but not specified'.format(key))
    sig, ra, dec, mod, r, n_obs, n_model = [params[key] for key in keys]
    field_density = round(n_model/(np.pi * (r*60)**2), 4)

    region = simple.survey.Region(survey, ra, dec) 
    print('Loading data...')
    stars = region.get_data('stars')
    galaxies = region.get_data('galaxies')
    print('Found {} stars...'.format(len(stars)))
    print('Found {} galaxies...'.format(len(galaxies)))
    iso = region.survey.get_isochrone(params['mod'])
    g_radius = get_g_radius(region, stars, iso)

    print('Making diagnostic plots for (SIG, {}, {}) = ({}, {}, {})...'.format(survey.catalog['basis_1'], survey.catalog['basis_2'], sig, ra, dec))
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    density_plot(axs[0][0], region, stars, g_radius, iso, 'stars')
    density_plot(axs[1][0], region, stars, g_radius, iso, 'blue_stars')
    density_plot(axs[2][0], region, galaxies, g_radius, iso, 'galaxies')

    star_plot(axs[0][1], region, stars, iso)
    star_plot_aperture(axs[0][2], region, stars, iso, r)

    cm_plot(axs[1][1], region, stars, iso, g_radius, 'stars')
    hess_plot(axs[1][2], region, stars, iso, g_radius)
    
    cm_plot(axs[2][1], region, galaxies, iso, g_radius, 'galaxies')
    radial_plot(axs[2][2], region, stars, galaxies, iso, g_radius, field_density)

    """
    # Name
    try: # ugali
        association_string = candidate['NAME']
    except: # simple
        # Check for possible associations
        glon_peak, glat_peak = ugali.utils.projector.celToGal(ra, dec)
        catalog_array = ['McConnachie15', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'WEBDA14','ExtraDwarfs','ExtraClusters']
        catalog = ugali.candidate.associate.SourceCatalog()
        for catalog_name in catalog_array:
            catalog += ugali.candidate.associate.catalogFactory(catalog_name)

        idx1, idx2, sep = catalog.match(glon_peak, glat_peak, tol=0.5, nnearest=1)
        match = catalog[idx2]
        if len(match) > 0:
            association_string = '{} at {:0.3f} deg'.format(match[0]['name'], float(sep))
        else:
            association_string = 'No association within 0.5 deg'

    association_string = str(np.char.strip(association_string))
    association_string = repr(association_string)
    """
    info_string = r'($\alpha$, $\delta$, $\mu$) = ({:0.2f}, {:0.2f}, {:0.2f})'.format(ra, dec, mod)
    detect_string = r'($\sigma$, $r$, n_obs, n_model) = ({:0.2f}, {:0.2f}, {:0.2f}, {:0.2f})'.format(sig, r*60, n_obs, n_model)

    #plt.suptitle(association_string+'\n' + info_string+'\n' + detect_string, fontsize=24)
    plt.suptitle(info_string+'\n' + detect_string, fontsize=24)

    file_name = 'candidate_{:0.2f}_{:0.2f}_{:0.2f}'.format(sig, ra, dec)
    plt.savefig(survey.output['plot_dir']+'/'+file_name+'.png',  bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True, help='config file')
    parser.add_argument('--sig_cut', type=float, required=False, help='Significance threshold to plot', default=5.5)
    parser.add_argument('--infile', type=str, required=False, help='Candidate list to plot. No further arguments necessary if passed')
    parser.add_argument('--sig', type=float, required=False, help='Candidate significance')
    parser.add_argument('--ra', type=float, required=False, help='RA')
    parser.add_argument('--dec', type=float, required=False, help='Dec')
    parser.add_argument('--r', type=float, required=False, help='Aperture')
    parser.add_argument('--modulus', type=float, required=False, help='Distance modulus')
    parser.add_argument('--n_obs', type=float, required=False, help='Observed stars')
    parser.add_argument('--n_model', type=float, required=False, help='Expected number of field stars')
    args = vars(parser.parse_args())

    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        survey = simple.survey.Survey(cfg)

    if args['infile']:
        try:
            candidate_list = np.load(args['infile'])
        except IOError:
            try:
                candidate_list = fitsio.read(args['infile'])
            except IOError:
                raise IOError('Infile not found or format not supported. Needs to be .npy or .fits')

        try: # simple
            candidate_list = candidate_list[candidate_list['SIG'] > args['sig_cut']]
        except: # ugali
            candidate_list = candidate_list[candidate_list['TS'] > args['sig_cut']]

        for candidate in candidate_list:
            make_plot(survey, candidate)

    else:
        make_plot(survey, sig=args['sig'], ra=args['ra'], dec=args['dec'], r=args['r'], mod=args['modulus'], n_obs=args['n_obs'], n_model=args['n_model'])
         
