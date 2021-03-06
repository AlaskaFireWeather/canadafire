import numpy as np
import canadafire

#lons = np.array([-147.7164])
lons = np.array([0., -147.7164])

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
msn,asn = canadafire.solar_noon(2019,-1,-1,lons,doy=1)    # should be same as jan 1
vals.append((2019,1,1,lons,msn,asn))
for Y in range(2019,2020):
    for D in range(1,2):
        for M in range(1,13):
            msn,asn = canadafire.solar_noon(Y,M,D,lons)
            vals.append((Y,M,D,lons,msn,asn))
#            print('{}-{}-{}: {} {}'.format(Y,M,D,msn,asn-msn))

for Y,M,D,lons,msn,asn in vals:
    diff_h = asn-msn
    diff_s = 3600*diff_h
    min = diff_s // 60
    sec = diff_s - min*60
    print('{}-{}-{}: {} {} {}:{}'.format(Y,M,D,msn,asn,min,sec))
