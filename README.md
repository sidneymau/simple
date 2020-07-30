# Simple Binning Suite 

The simple binning suite includes modules for peforming a simple binning search of stellar overdensities in DES &c data, compiling a candidate list from the results of the search, and producing diagnostic plots from a candidate list.

This code has been adapted from [Keith Bechtol](https://github.com/bechtol)'s original simple binning program to be more modular and ready-to-use for future searches.

Maintained by [Sidney Mau](https://github.com/SidneyMau) and [Mitch McNanna](https://github.com/mcnanna). 

## Dependencies

This program uses the following:

python;
numpy,
scipy,
healpy,
astropy,
matplotlib, 
[fitsio](https://github.com/esheldon/fitsio),
[pyyaml](https://pyyaml.org/)

[ugali](https://github.com/DarkEnergySurvey/ugali),
[batchtools](https://github.com/kadrlica/batchtools) (for farming).

## Configuration and use

`simple` expects the catalog data to be in a specific format: each data file should be a .fits file corresponding to a single healpix pixel with `nest=False`. 

You will want to create a directory such as `simple_run/` in which to perform the search. 

Running `init_simple.py` in this directory will produce output directories as well as a `config.yaml` file.
This `config.yaml` file specifies the configuration for the search and should be edited for your use case. The details of this config file are given their own section below.

Running `search.py --config config.yaml --ra RA --dec DEC --outfile results.npy` will perform the search in the healpix pixel containing the specified RA and DEC, writing the results to `results.npy` file. 
To run a parallel search over the entire dataset using the HTCondor batch system, use `parallel_search.py --config config.yaml`. Output will be written to the `results_dir` specified in the config file. 
The results can then be compiled into a single candidate list FITS file by running `make_list.py`

Diagnostic plots of a list of candidates can be made using `plotting/diagnostic_plots.py --config config.yaml --infile results.npy`. Unless the optional `--sig_cut` argument is passed, only candidates with `SIG` > 5.5 will be plotted. 
`plotting/farm_plots.py --config config.yaml --infile results.npy` can be used to parallelize the plotting using HTCondor. Plots will be written to the `plot_dir` specified in the config file.  

## Config File

The fields of the config file are explained below. *Italicised* field names can be set to `null` if not applicable. *(Parenthetical)* field names can be entirely left out of the config file. 

* `band_1`: Isochrone selection will be done in `band_1 - band_2` vs `band_1` color space. Ex: `G` 
* `band_2`: Isochrone selection will be done in `band_1 - band_2` vs `band_1` color space. Ex: `R` 
* *`band_3`*: An additional selection will be done in `band_2 - band_3` vs `band_2` color space. Ex: `I` 
* `survey`:   
  * *`fracdet`*: Coverage area fraction map, used in background density calculations  
* `catalog`: 
  * `dirname`: Path to directorty containing healpixelized data files 
  * `nside`: nside of healpixel 
  * `mag_max`: Magnitude limit for stars included in the search
  * `basis_1`: Longitudinal coordinate, e.g. `RA`
  * `basis_2`: Latitudinal coordinate, r.g. `DEC`
  * `mag`: Name of magnitude field in the .fits file, with the band replaced by `{}`. Ex: `PSF_MAG_{}`
  * `mag_err`: Name of magnituded error field in the .fits file, with the band replaced by `{}`. EX: `PSF_MAG_ERR_{}`
  * *`reddening`*: Reddening correction to be subtracted from `mag`, with the band replaced by `{}`: `mag_deredenned = mag - reddening`. Ex: `EBV_{}`. 
  * *`quality`*: Quality selection to reduce data. Ex: `PSF_MAG_G < 24.5 && PSF_MAG_R < 24.0`
  * *`stars`*: Selection for stars. Ex: `EXTENDED_CLASS < 1` (Note: required for plotting)
  * *`galaxies`*: Selection for galaxies. Ex:`EXTENDED_CLASS >= 1` (Note: required for plotting)
  * *`other`*: An additional selection contained in a `.py` file. The file must contain a function named `sel(survey, data)` which takes a `survey.Survey` object and data array as arguments, and returns a boolean array of the same length as `data`. Ex: `extra_selection`, where there is a file named `extra_selection.py` in the directory where the search is being run. 
  * *`size`*: Name of morphological size parameter
* `isochrone`:
  * `name`: Name of `ugali` isochrone class to use
  * `survey`: Name of survey to use in isochrone creation
  * `age`: Isochrone age in Gyr 
  * `metallicity`: Isochrone metallicity Z. 
* `output`:
  * `results_dir`: Directory to store search results files
  * `log_dir`: Directory to store log files
  * `plot_dir`: Directory to store plots
* `batch`:
  * `max_jobs`: Maximum jobs to run simultaneously on HTCondor

The following fields are entirely optional and can all be left out.
* *`(grid)`*:
  * *`(delta_x)`*: Step size for star binning in degrees. Default: 0.01
  * *`(smoothing)`* Gaussian smoothing kernal size in arcmin. Default: 2.0
  * *`(bin_edge)`* (Square) grid radius in degrees. Default: 8.0
  * *`(grid_dir)`* Directory containing pre-made grid array selections.
* *`(moduli)`*:
  * *`(start)`*: Closest distance modulus to search. Default: 16
  * *`(stop)`*: Farthest distance modulus to search. Default: `mag_max`
  * *`(step)`*: Step size between distance moduli. Default: 0.5 
 

## Notes

* By default, the plotting scripts will only produce plots for hotspots with statistical significance greater than 5.5 sigma. This threshold has intentionally chosen to be low (such that investigations can be made of very low-significance hotsposts and candidates) but also to minimize junk. This can be changed by specifying the optional `--sig_cut` argument to the plotting scripts.


* One of the most computationally intensive steps in the search algorithm is the construction of a 2D boolean array representing which coordinates on a fine 2D grid are included in a given healpixel. This becomes a significant bottlneck when the `delta_x` parameter in the `grid` config is lowered below the default `0.01`. In such a case, the creation of these arrays can be done as a pre-processing step by using `create_grid_selections.py`. 
