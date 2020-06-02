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

Running `search.py --config config.yaml --ra RA --dec DEC --outfile out.npy` will perform the search in the healpix pixel containing the specified RA and DEC, writing the results to a `.npy` file. To run a parallel search over the entire dataset using the HTCondor batch system, use `parallel_search --config config.yaml`. 





`config.yaml` handles most of the setup and can be modified according to use.

You will need to create a directory such as `simple_run/`.
Copy `config.yaml` into `simple_run/`; this will let the the scripts know where to source the config info from as well as where to write the output.
If you have fracdet or maglim files, those should also be copied or linked here; specify their filepaths in `config.yaml`.

To run the simple binning search, run `farm_simple.py` in `simple_run/`.
This will run `search_algorithm.py` over the given data set and write the output to `results_dir/`, logging each job in `results_dir/log_dir`.
The results can then be compiled into a candidate list by running `make_list.py` from `simple_run/` (this is saved as a `.npy` file).

To produce plots from a candidate list, run `farm_plots.py` in `simple_run/`.
The output will be written to `save_dir/` with logs in `save_dir/log_dir/`.

## Notes

By default, `farm_plots.py` will only produce plots for hotspots with statistical significance greater than 5.5 sigma.
This threshold has intentionally chosen to be low (such that investigations can be made of very low-significance hotsposts and candidates) but also to minimize junk.
