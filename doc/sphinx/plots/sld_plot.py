import pylab
import periodictable.nsf

params = {'figure.figsize': (10,8)}
pylab.rcParams.update(params)

periodictable.nsf.sld_plot()
