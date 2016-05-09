# P3R
P3R stands for: **Python PhreeqC Parameter Refinement**. It's a Python-PhreeqC coupled code for optimization and statistical analyses of chemical parameters. 

Only Author so far is Frank Heberling, a Geochemist from Karlsruhe. Any contributions / ideas / suggestions are welcome. So if you want to join in, just leave me a note.

The GUI and the optimization routine are largely adapted from the sxrd modul in the **xraypy** project on https://github.com/xraypy/tdl, even though P3R contains quite some new features.

The chemical equilibrium calculations fully rely on **PhreeqC** (http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/). Only Parameter variations, optimization and statistical analyses are performed by phyton code. Correspondingly "Chemical Parameters" from the title means any type of value, which can be input into a PhreeqC calculation.

PhreeqC is implemented into the code through **PhreeqPy** (http://www.phreeqpy.com/).

Instalation/Dependencies:

To run the code you need: 
    **Python** (https://www.python.org/,  only version 2.7 is used so far, but in principle it should run on other versions as well).
    **Numpy** (http://www.numpy.org/)
    **Scipy** (https://www.scipy.org/)
    **matplotlib** (http://matplotlib.org/)
    **wxpython** (http://www.wxpython.org/)
    in order to install PhreeqPy I recommend to get **pip** (https://bootstrap.pypa.io/get-pip.py; then you do python get-pip.py; after you have installed pip you can run: pip install -U phreeqpy...; copy a suitable IPhreeqc.dll from here: ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/; then you should be ready to run P3R)
    
    Further instructions and some example P3R projects will follow...


