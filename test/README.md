

Main Reference in Van Wagner et al, 1985:
    Canadian Forest Service. Forestry Technical Report
    https://cfs.nrcan.gc.ca/publications?id=19973
    0-662-13906-2

VanWagner1985SampleOfOutput.csv:
    Published results including test input and output, digitized from
    Van Wager 1985 (p 17)

canadafire_w85.f:
    Published FORTRAN 77 code, digitzed from Van Wagner 1985 (p 11-15)


f32in.dat:
    Published test input file, digitized from Van Wagner 1985 (p 16)

f32out.dat: (GENERATED FILE)
    Result of running canadafire_w85.f on f32in.dat.
    SHOULD MATCH VanWagner1985SampleOfOutput.csv

canadafire.f:
    Fortran code received by Elizabeth Fischer in August 2021
    (Python extension C code is a translation of this).

canadafire_2021.f:
    Driver code to compare canadafire.f directly with canadafire_w85.f
