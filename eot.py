import numpy as np
import canadafire

lons = np.array([-147.7164])
#lons = np.array([0.])

dates = [
    (2019,1,1),
    (2020,1,1),
]
#    (2020,4,1),
#    (2020,7,1),
#    (2020,9,1),
#    (2020,10,1),
#]

vals = list()
for Y in range(2019,2020):
    for D in range(1,2):
        for M in range(1,13):
            msn,asn = canadafire.solar_noon(Y,M,D,lons)
            vals.append((Y,M,D,lons,msn,asn))
#            print('{}-{}-{}: {} {}'.format(Y,M,D,msn,asn-msn))

for Y,M,D,lons,msn,asn in vals:
    print('{}-{}-{}: {} {}'.format(Y,M,D,msn,asn))
