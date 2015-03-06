# Radio Weak Lensing Simulations

Creates an empty measurement set (CASA Table; MS) then simulates visibilites into it, given a sky model.
The MS paremeters can be set via in the `pyxis-rwlsims.conf` file.

This repo comes with a MeerKAT antenna table (MeerKAT64_ANTENNAS), a MeqTrees TDL profiles file (tdlconf.profiles) and a sky model (gauss.txt). 

The main functions are `make_empty_ms` and `simsky`. Both can be called from the command line as follows:  

```
$ pyxis msname function[parameter=value]
```
Run `$ pyxis help[function]` for help on any of these functions.

## Example
```
$ pyxis meerkat_trial1.MS OUTDIR=Test LSM=gauss.txt make_empty_ms simsky im.lwimager.make_image[restore=True,operation=csclean]
```

The above function will create an empty MS ,*meerkat_trial1.MS*, simulate the sky model *gauss.txt* into it then make a clean map using  **lwimager** (if its installed in your system). The images will placed in *OUTDIR/plots-foo/*. 

