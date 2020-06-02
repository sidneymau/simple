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

`simple` expects the catalog data to be in a specific format: each data file should correspond to a single healpix pixel with `nest=False`. 

You will want to create a directory such as `simple_run/` in which to perform the search. 

Running `init_simple.py` in this directory will produce output directories as well as a `config.yaml` file.
This `config.yaml` file specifies the configuration for the search and should be edited for your use case. The details of this config file are given their own section below.

Running `search.py --config config.yaml --ra RA --dec DEC --outfile results.npy` will perform the search in the healpix pixel containing the specified RA and DEC, writing the results to `results.npy` file. 
To run a parallel search over the entire dataset using the HTCondor batch system, use `parallel_search.py --config config.yaml`. Output will be written to the `results_dir` specified in the config file. 
The results can then be compiled into a single candidate list FITS file by running `make_list.py`

Diagnostic plots of a list of candidates can be made using `plotting/diagnostic_plots.py --config config.yaml --infile results.npy`. Unless the optional `--sig_cut` argument is passed, only candidates with `SIG` > 5.5 will be plotted. 
`plotting/farm_plots.py --config config.yamls --infile results.npy` can be used to parallelize the plotting using HTCondor. Plots will be written to the `plot_dir` specified in the config file.  


## Notes

By default, the plotting scripts will only produce plots for hotspots with statistical significance greater than 5.5 sigma. This threshold has intentionally chosen to be low (such that investigations can be made of very low-significance hotsposts and candidates) but also to minimize junk. This can be changed by specifying the optional `--sig_cut` argument to the plotting scripts.
