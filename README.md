# halphagui
development of gui interface for H-alpha narrowband imaging


# some issues I ran into setting up ginga

when installing ginga using pip (anaconda 3.7), the util module is out of date.

* I downloaded the ginga source code from github, and then ran setup.py
* I then copied the util directory to my anaconda installation of ginga
  * cp -r ~/github/ginga/ginga/util ~/anaconda3/lib/python3.7/site-packages/ginga/.
* it's probably better to copy the entire package, but this is what I did for starters.

# some issue with getting qt to run

at first, I wasn't able to run any of the ginga examples, and this was because I hadn't defined an environment variable correctly.  Sorry - I can't remember the details now, but the instructions set to set the environment variable to qt, but it actually needed to be qt5.  

