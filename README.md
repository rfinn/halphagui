# halphagui
- development of gui interface for the analysis of H-alpha narrowband imaging of nearby galaxies
- the gui was developed using HDI imaging data from the KPNO 0.9m.  Will also test on KPNO Mosaic data and INT WFC data.

# Workflow #

Edit the gui using qt designer:
```
designer &
```
Then open a particular ui file.  
```
pyuic5 -x halpha.v1.ui -o halphav1.py
python halphamain.py
```

# Requirements ##
see wiki for more detailed instructions
## Python Modules (need to add versions) ##


* ginga
* astropy
* pyqt5
* photutils
* statmorph https://statmorph.readthedocs.io/en/latest/

