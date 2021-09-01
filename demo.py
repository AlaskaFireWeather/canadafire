import numpy as np
import inspect
import canadafire

# Sample program using the canadafire Python extension.

print(canadafire.canadafire.__doc__)
print()

tin = np.array( [[15.,16.], [31.,32.]] )
hin = np.array( [[70.,80.], [20.,15.]] )
win = np.array( [[10.,5.], [60.,80.]] )
rin = np.array( [[30.,40.], [0., 1.]] )
imonth = np.array([7,8], dtype='i')

bui,ffmc,isi,fwi,dmc,dc = canadafire.canadafire(tin,hin,win,rin,imonth, debug=True)

print('\nbui:')
print(bui)

print('\nffmc:')
print(ffmc)

print('\nisi:')
print(isi)

print('\nfwi:')
print(fwi)

print('\ndmc:')
print(dmc)

print('\ndc:')
print(dc)
