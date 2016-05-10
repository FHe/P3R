# P3R
**P3R** stands for: **Python PhreeqC Parameter Refinement**. It's a Python-PhreeqC coupled code for optimization and statistical analyses of chemical parameters. 

Author is Frank Heberling. Any contributions / ideas / suggestions are welcome. So if you want to join in developing, just leave me a note.

The GUI and the optimization routine are largely adapted from the **sxrd modul** in the **TDL** project on https://github.com/xraypy/tdl, although P3R contains some stunning new features.

The chemical equilibrium (and other) calculations fully rely on **PhreeqC** (http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/). Correspondingly "chemical parameters" in the title means any type of value, which can be input into a PhreeqC calculation. Visualization, Parameter variations, optimization and statistical analyses are performed by phyton code.

PhreeqC is implemented in P3R through **PhreeqPy** (http://www.phreeqpy.com/).

***Installation / Dependencies:***

To run the code you need: 
   *  **Python** (https://www.python.org/,  version 2.7 is used so far, but there is no obvious reason why it should not run on other versions).
   *  **Numpy** (http://www.numpy.org/)
   *  **Scipy** (https://www.scipy.org/)
   *  **matplotlib** (http://matplotlib.org/)
   *  **wxpython** (http://www.wxpython.org/)
   * in order to install **PhreeqPy** I recommend to get **pip** (https://bootstrap.pypa.io/get-pip.py and execute **python get-pip.py**; (if you don't have it already with your Python installation); after you have installed pip you can run: **pip install -U phreeqpy**; and copy a suitable **IPhreeqc.dll** from here: ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/) 
   * now you should be ready to run **P3R**
    
Further instructions and some example P3R projects will follow...


