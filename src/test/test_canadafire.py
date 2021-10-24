import os
import numpy as np
import canadafire
import pytest

"""Unit test for canadafire C module.  To run test:
     pytest test_canadafire.py
"""

def read_data(fname):
    """Read raw data files of inputs and outputs,

    Structure should be [number of ensembles x 183]
    For the year 2018, ensemble number = 51

    NOTE: The 183 columns are split between 4 lines of 50 columns
          each; + 33 on the last line of 4.

    Returns: np.array[num_ensembles, 183 days]
    """
    data = list()
    with open(fname, 'r') as fin:
        while True:
            try:
                sline = \
                    next(fin).split() + \
                    next(fin).split() + \
                    next(fin).split() + \
                    next(fin).split()
                data.append([float(x) for x in sline])
            except StopIteration:
                break

    ret = np.array(data)
    return ret

def compare(vname, val0, val1, tolerance):
    """Find number of differences between two Numpy arrays"""
    diff = val1 - val0
    ndiff = np.sum(np.abs(diff) >= tolerance)
    print('{}: {}'.format(vname, ndiff))
    return ndiff


def test_canadafire():

    dir = os.path.dirname(__file__)

    # Sample Inputs
    tin = read_data(os.path.join(dir, '2018_tmp2m_anom_wgt.txt'))
    hin = read_data(os.path.join(dir, '2018_rh_anom_wgt.txt'))
    win = read_data(os.path.join(dir, '2018_uv_wgt.txt'))
    rin = read_data(os.path.join(dir, '2018_precip_anom_wgt.txt'))

    # Previously computed outputs for sample inputs
    bui0 = read_data(os.path.join(dir, '2018_bui_anom_wgt.txt'))
    ffmc0 = read_data(os.path.join(dir, '2018_fm_anom_wgt.txt'))
    isi0 = read_data(os.path.join(dir, '2018_isi_anom_wgt.txt'))
    fwi0 = read_data(os.path.join(dir, '2018_fwi_anom_wgt.txt'))
    dmc0 = read_data(os.path.join(dir, '2018_dmc_anom_wgt.txt'))
    dc0 = read_data(os.path.join(dir, '2018_dc_anom_wgt.txt'))

    # Construct the standard imonth array (OPTIONAL)
    imonth = np.array([4]*30 + [5]*31 + [6]*30 + [7]*31 + [8]*31 + [9]*30, dtype='i')

    bui1,ffmc1,isi1,fwi1,dmc1,dc1 = canadafire.canadafire(tin,hin,win,rin)

    # Compare previously computed output to what we got now.
    # Files have 3 decimal places of precision
    assert 0 == compare('bui', bui0, bui1, 1.e-3)
    assert 0 == compare('ffmc', ffmc0, ffmc1, 1.e-3)
    assert 0 == compare('isi', isi0, isi1, 1.e-3)
    assert 0 == compare('fwi', fwi0, fwi1, 1.e-3)
    assert 0 == compare('dmc', dmc0, dmc1, 1.e-3)
    assert 0 == compare('dc', dc0, dc1, 1.e-3)


    # ----------------------------------------
    # Test the typechecking
    with pytest.raises(TypeError):
        bui1,ffmc1,isi1,fwi1,dmc1,dc1 = canadafire.canadafire(tin.astype('f'),hin,win,rin)
