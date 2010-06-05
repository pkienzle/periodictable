import pylab
import periodictable
fig_width = 10
fig_height = 8
fig_size =  [fig_width,fig_height]
params = {'figure.figsize': fig_size}
pylab.rcParams.update(params)
Fe_2 = periodictable.Fe.ion[2]
Q = pylab.linspace(0,16,200)
M = Fe_2.magnetic_ff[Fe_2.charge].j0_Q(Q)
pylab.xlabel(r'Magnetic Form Factor for Fe')
pylab.ylabel(r'$\AA^{-1}$')
pylab.title('Ion specific property for Fe')
pylab.plot(Q,M)
