import periodictable
import pylab

params = {'figure.figsize': (10,8)}
pylab.rcParams.update(params)

D = [(el.number,el.density,el.symbol)
    for el in periodictable.elements]

bbox = dict(boxstyle="round",lw=1,ec=(0,0,0),fc=(0.85,0.8,0.8))
for Z,density,sym in D:
    if density is not None:
        pylab.text(Z,density,sym,bbox=bbox)
pylab.axis([0,110,0,25])
pylab.xlabel('Element number')
pylab.ylabel('Density of element')
pylab.title('Density for elements')
