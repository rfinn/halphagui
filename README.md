
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

# Setting up a virtual environment
Check that you are running some version of python3.

To create a virtual environment,
```
cd ~/github/halphagui
```

Then create the environment.  We are calling it venv
```
python -m venv venv
```

To activate the new environment
```
source venv/bin/activate
```
To install the required modules, type:

```
pip install -r requirements.txt
```
To deactivate the virtual environment
```
(venv) $ deactivate
```

# Requirements ##

see wiki for more detailed instructions


## Python Modules (need to add versions) ##


* ginga
* astropy
* pyqt5
* photutils
* statmorph https://statmorph.readthedocs.io/en/latest/

