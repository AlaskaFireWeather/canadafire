import numpy as np
import canadafire

#lons = np.array([-147.7164])
lons = np.array([0.])

dates = [
    (2019,1,1),
    (2020,1,1),
]
#    (2020,4,1),
#    (2020,7,1),
#    (2020,9,1),
#    (2020,10,1),
#]

for Y in range(2019,2020):
    for D in range(1,3):
        for M in range(1,2):
            msn,asn = canadafire.solar_noon(Y,M,D,lons)

            print('{}-{}-{}: {} {}'.format(Y,M,D,msn,asn-msn))
