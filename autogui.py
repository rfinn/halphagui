#!/usr/bin/env python

import os

# 2017 data
pointings = ['1','3','4','6','7','10','12','13','14','16','19','23','27','28']
prefix = 'v17'
obsyear = '2017'
for p in pointings:
    pprefix = prefix+'p%02d'%(int(p))
    os.system('python ~/github/halphagui/testing/halphamain.py --virgo --laptop --pointing '+p+' --auto --prefix '+pprefix+' --obsyear '+obsyear)
