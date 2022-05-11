import subprocess
import os
import re
import pandas as pd
import canadafire
import functools
import contextlib

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

header_re = re.compile( (r'\s+([^\d\s]+)'*12) + r'\s*')
data_re = re.compile( (r'\s+([\d\.]+)' * 13) + r'\s*')


# -------------------------------------------
def read_f32out(fname):
    """Reads the output of the Fortran test programs"""

    # Parse the output file
    rows = list()
    with open(fname,'r') as fin:
        for line in fin:
            line = line.rstrip()

            # See if it's a header line
            match = header_re.match(line)
            if match is not None:
                # Replace "DATE" with "MONTH" and "DAY"
                headers = ['MONTH', 'DAY'] + [match.group(i) for i in range(2,13)]
                continue

            # See if it's a data row
            match = data_re.match(line)
            if match is not None:
                vals = [int(match.group(1)), int(match.group(2))] + \
                    [float(match.group(i)) for i in range(3,14)]

                # Add to output
                rows.append(dict(zip(headers, vals)))
                continue

    # Assemble dataframe of what the legacy FORTRAN 77 code produced.
    df1 = pd.DataFrame(rows)
    return df1
# -------------------------------------------

def get_fortran(suffix):
    try:
        os.chdir(TEST_DIR)
        cmd = ['gfortran', f'canadafire{suffix}.f', '-o', f'canadafire{suffix}']
        subprocess.run(cmd, check=True)

        cmd = [os.path.join(TEST_DIR, f'canadafire{suffix}')]
        subprocess.run(cmd, check=True)

        df1 = read_f32out(f'f32out{suffix}.dat')

        return df1

    finally:
        with contextlib.suppress(FileNotFoundError):
            os.remove(f'f32out{suffix}.dat')
            os.remove(f'canadafire{suffix}')


# ------------------------------------
def _col(df, cname):
    val = df[cname].to_numpy().astype('d')
    val = val.reshape(1,len(val))
    return val

@functools.lru_cache()
def get_answer():
    df0 = pd.read_csv('VanWagner1985SampleOfOutput.csv')
    #df0 = df0.drop('DSR', axis=1)
    return df0

@functools.lru_cache()
def get_pythonext():
    """Tests the C-based Python extension"""

    df0 = get_answer()
    tin = _col(df0, 'TEMP')
    hin = _col(df0, 'RH')
    win = _col(df0, 'WIND')
    rin = _col(df0, 'RAIN')
    imonth = df0.MONTH.to_numpy().astype('i')

    ocf = canadafire.canadafire(
        tin, hin, win, rin,
        imonth=imonth,
        ffmc0=85.5, dmc0=6.0, dc0=15.0)

    ucols0 = ['BUI', 'FFMC', 'ISI', 'FWI', 'DMC', 'DC']
    ucols1 = ['FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI']
    outs0 = {colname:x[0,:] for colname,x in zip(ucols0,ocf)}

    # Reorder columns to standard order for comparison
    outs = {'MONTH':imonth, 'DAY':0, 'TEMP':tin.reshape(-1), 'RH':hin.reshape(-1), 'WIND':win.reshape(-1), 'RAIN':rin.reshape(-1)}
    for x in ucols1:
        outs[x] = outs0[x]
    df2 = pd.DataFrame(outs)
    return df2

def compare_answer(df1, df0, tolerance):
    # Read what the answer SHOULD be
#    df0 = pd.read_csv('VanWagner1985SampleOfOutput.csv')

    # Assert that ALL values in the two dataframes are equal
    for cname in ('BUI', 'FFMC', 'ISI', 'FWI', 'DMC', 'DC'):
        # Only avoid if BOTH don't have this column.
        if cname not in df0 and cname not in df1:
            continue
        assert ((df1[cname] - df0[cname]).abs() < tolerance[cname]).all(), \
            f"Columns {cname} don't match"

# ==================================================================

# Original Fortran code must match published answer exactly
tolerance_w85 = dict(FFMC=.001, DMC=.001, DC=.001, ISI=.001, BUI=.001, FWI=.001, DSR=.001)
# Derivative Fortran code can be a little off
tolerance_2021 = dict(FFMC=0.1, DMC=0.001, DC=0.001, ISI=0.3, BUI=0.001, FWI=0.4, DSR=0.1)

# C code should match Derivative Fortran code closely
tolerance_pyext = dict(FFMC=0.06, DMC=0.06, DC=0.06, ISI=0.06, BUI=0.06, FWI=0.06, DSR=0.06)
tolerance_pyext_df0 = dict(FFMC=0.1, DMC=0.06, DC=0.06, ISI=0.3, BUI=0.06, FWI=0.4, DSR=0.1)

def test_fortran_w85():
    print('------------------------ test_fortran_w85')
    compare_answer(get_fortran('_w85'), get_answer(), tolerance_w85)

def test_fortran_2021():
    print('------------------------ test_fortran_2021')
    df0 = get_answer()
    df1 = get_fortran('_2021')
    compare_answer(df1, df0, tolerance_2021)


def test_pythonext():
    df1 = get_fortran('_2021').drop('DSR', axis=1)
    df0 = get_answer().drop('DSR', axis=1)
    df2 = get_pythonext()

    # C code should match Derivative Fortran code closely
    compare_answer(df2, df1, tolerance_pyext)

    # ...and should match original Van Wagner 1985 answer with same
    # tolerance as the derivative Fortran code.
    compare_answer(df2, df0, tolerance_pyext_df0)

