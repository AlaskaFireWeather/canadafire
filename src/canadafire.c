#define PY_SSIZE_T_CLEAN
#include <Python.h>                // Must be first
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


static char module_docstring[] = "\
Module computes the Canada Fire Index on a gridded dataset.";

// For help with Python extension:
// https://stackoverflow.com/questions/56182259/how-does-one-acces-numpy-multidimensionnal-array-in-c-extensions



// =====================================================================
// Helper Functions

static inline double sqr(double x) { return x*x; }

/* computes base^8 */
static inline double ipow8(double base)
{
    double const base2 = base*base;
    double const base4 = base2*base2;
    double const base8 = base4*base4;
    return base8;
}

/** Computes exponentiation using repeated squares: base^exp */
static inline double ipow(double base, int exp)
{
    double result = 1;
    for (;;) {
        if (exp & 1) result *= base;
        exp >>= 1;
        if (!exp) break;
        base *= base;
    }

    return result;
}

/** Checks that arr has same rank and dimensions as arr0 */
bool check_rank_dims(char const *name0, PyArrayObject *arr0, char const *name, PyArrayObject *arr)
{
    char msg[256];

    // Check rank
    if (PyArray_NDIM(arr) != PyArray_NDIM(arr0)) {
        sprintf(msg, "Parameter %s must have same rank as %s", name, name0);
        PyErr_SetString(PyExc_TypeError, msg);
        return false;
    }

    // Check dimensions
    for (int i=0; i<PyArray_NDIM(arr); ++i) {
        if (PyArray_DIM(arr,i) != PyArray_DIM(arr0,i)) {
            sprintf(msg, "Parameter %s must have same dimensions as %s: but dimension %d is %ld instead of %ld", name, name0, i, (long)PyArray_DIM(arr,i), (long)PyArray_DIM(arr0,i));
            PyErr_SetString(PyExc_TypeError, msg);
            return false;
        }
    }

    return true;
}


// =====================================================================
// The core computation, converted from Fortran

// this is twice the value in instructions
static const double le[] = {6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0};
static const double lf[] = {-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6};


static void canadafire(
double const t, double const h, double const w, double const r,    // Input values
int im, // The current month
double const ffmc0, double const dmc0, double dc0,    // Initial values (accumulated from yesterday)
double *_bui, double *_ffm, double *_isi, double *_fwi, double *_dsr, double *_dmc, double *_dc)
{
    bool iswitch;

    // EQ 1: fine fuel moisture code 
    double mo = (147.2*(101.-ffmc0))/(59.5+ffmc0);

    // EQ 2: calculate rf
    double rf = r;
    if (rf > 0.5) {
        rf = r - 0.5;

        // If r<0.5 then you can not use the following equations 3a and 3b
        // then the rainfall routine must be omitted
        // Condition 1 
        // calculate wmo (fine fuel moisture code from previous day, mo) 
        double mr;
        if (mo <= 150.0) {
            // EQ 3a
            mr = mo
                + 42.5*rf*exp(-100./(251.-mo)) * (1.-exp(-6.93/rf));
        } else {
            // EQ 3b
            mr = mo
                + 42.5*rf*exp(-100./(251.-mo)) * (1.-exp(-6.93/rf))
                + 0.0015 * sqr(mo-150.) * sqrt(rf);
        }

        // Condition 3
        if (mr > 250.0) mr = 250.0;

        mo = mr;
    }

    // EQ 4: calculate ed
    double ed = 0.942 * pow(h,0.679)
        + 11. * exp((h-100.)/10.)
        + 0.18 * (21.1-t) * (1.-1./exp(0.115*h));

    // iswitch = false;

    double m;
    if(mo > ed) {
        // EQ 6a
        double const ko = 0.424*(1. - pow((h/100.0), 1.7))
            + (0.0694 * sqrt(w))* (1. - ipow8(h/100.));

        // EQ 6b
        double const kd = ko* 0.581*exp(0.0365*t);

        // EQ 8
        m = ed+(mo-ed) * pow(10.0, (-1.*kd));
        // iswitch = true;
    } else {
        // EQ 5 
        double const ew = 0.618 * pow(h, 0.753)
            + 10.*exp((h-100.)/10.)+0.18*(21.1-t) * (1.-exp(-0.115*h));

        if(mo < ew) {
            // EQ 7a
            double const kl =
                0.424 * (1.- pow((100.-h)/100, 1.7))
                + (0.0694 * sqrt(w)) * (1. - ipow8((100.-h)/100.));

            // EQ 7b
            double const kw = kl * (0.581 * exp(0.0365*t));

            // EQ 9 
            m = ew - (ew-mo)*pow(10.0, (-1.*kw));
        } else {
            m = mo;
        }
        iswitch = true;
    }

    // direction 8
    // if (!iswitch) m = mo;
    // Condition 2 
    // if(m > 250.) m = 250.;

    // calculate fine fuel moisture code
    // EQ 10
    double ffm = (59.5 * (250.-m)) / (147.2+m);

    // restrictions on ffm
    if(ffm> 101.) ffm = 101.;
    if(ffm <= 0.) ffm = 0.;

    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    // Duff moisture code
    //
    // yesterdays code becomes dmc0
    // direction 1
    double const po = dmc0;

    // condition 1 r must be above zero to calculate these quantities. 
    double pr;
    if (r <= 1.5) {
        pr = po;
    } else {
        // EQ 11
        double const re = 0.92*r - 1.27;

        // EQ 12
        mo = 20.0 + 280.0 / exp(0.023*po);

        //----------------------------------------------------------
        // calculate b depending on po value
        double b;
        if(po <= 33.) {
            // EQ 13a
            b = 100./(0.5+0.3*po);
        } else {
            // EQ 13c
            b = 6.2*log(po)-17.2;
            if(po-65.0 <= 0.) b = 14.-1.3*log(po);
        }

        // EQ 13b
        //----------------------------------------------------------
        // EQ 14
        double mr = mo + (1000.*re)/(48.77+b*re);

        // EQ 15
        pr = 43.43*(5.6348-log(mr-20.));
    }

    double k;
    if(t >= -1.1) {
        k=1.894*(t+1.1)*(100.-h)*(le[im-1]*0.0001);
    } else {
        k=0.0;
    }

    // condition 2 that pr can not be less than zero
    if (pr <= 0.) pr=0.;
    double dmc = pr+k;

    // condition 3 t must be greater than -1.1       
    if (dmc <= 0.) dmc=0.;
 
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    // drought code
    //
    // condition 1
    if(r > 2.8) {
        // EQ 18  
        double const rd = 0.83*r - 1.27;
        // EQ 19
        double const qo = 800.*exp(-dc0/400.);
        // EQ 20  
        double const qr = qo+3.937*rd;

        // EQ 21  
        double const dr = 400. * log(800./qr);
        if (dr > 0.0) {
           dc0 = dr;
        } else {
           dc0 = 0.0;
        }
    }

    // condition 3
    double v;
    if (t < -2.8) {
        v = lf[im-1];
    } else {
        // EQ 22
        v = 0.36*(t+2.8) + lf[im-1];
        // condition 4
    }
    if (v <= 0.) v=0.;

    // EQ 23
    double d = dc0 + 0.5*v;
    // set drought code
    double dc = d;

    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    // initial spread index, buildup index, fire weather index
    // initial spread index
    // EQ 24
    double const fw = exp(0.05039*w);
    // EQ 24
    double ff = 91.9*exp(-0.1386*m) * (1.+ pow(m,5.31)/4.93e07);
    // EQ 26
    // initial spread index
    double isi = 0.208 * fw * ff;

    // buildup index
    // EQ 27a
    double u;
    if(dmc <= 0.4*dc) {
        u = 0.8*dmc*dc/(dmc+0.4*dc);
    } else {
        // EQ 27b
        u = dmc -
            (1. -0.8*dc/(dmc+0.4*dc)) * (0.92 + pow(0.0114*dmc, 1.7));
    }

    // buildup index
    if (u < 0.) u = 0.0;

    // fire weather index
    // EQ 28a
    double fd;
    if(u <= 80.) {
        fd = 0.626*pow(u,0.809) + 2.;
    } else {
        // EQ 28b
        fd=1000./(25. + 108.64*exp(-0.023*u));
    }

    // EQ 29 
    double b= 0.1*isi*fd;

    // EQ 30a
    double s;
    if(b > 1.) {
        s = exp(2.72 * pow(0.434*log(b), 0.647));
    } else {
        s=b;
    }

    double const fwi = s;

    // daily severity rating
    double dsr = 0.0272 * pow(fwi,1.77);

    // ------------------------------------------
    // Return values
    *_bui = u;        // Buildup Index
    *_ffm = ffm;
    *_isi = isi;
    *_fwi = fwi;
    *_dsr = dsr;
    *_dmc = dmc;
    *_dc = dc;

}

// ===================================================================
// The Python-callable method

static char canadafire_canadafire_docstring[] = "\
Usage:\n\
    bui, ffmc, isi, fwi, dmc, dc = canadafire(tin, hin, win, rin, imonth, \n\
    ffmc0=85.5, dmc0=6.0, dc0=15.0, debug=False)\n\
\n\
Compute the canadafire index on a number of regions / gridcells.  See:\n\
    https://cwfis.cfs.nrcan.gc.ca/background/summary/fwi\n\
\n\
Dimensions\n\
----------\n\
  xy:\n\
    One or more space dimensions.  Typically x and y; or a\n\
    set of individual regions.  There is no interaction between\n\
    gridpoints.\n\
  time:\n\
    One data value per day.  The Canada Fire Index is defined for days\n\
    in April through September; although the calculation can be run\n\
    experimentally for other sets of days.\n\
\n\
NOTE: Optimized for row major (C order) arrays, with time dimension\n\
      having the smallest stride.\n\
\n\
Positional Arguments\n\
--------------------\n\
tin(xy,time): [degC] double\n\
    Temperature\n\
    (Picked at solar noon; or 2200UTC in Alaska)\n\
hin(xy,time): [%] double (0-100)\n\
    Relative humidity\n\
    (Picked at solar noon; or 2200UTC in Alaska)\n\
win(xy,time): [km/h] double\n\
    Wind speed\n\
    (Picked at solar noon; or 2200UTC in Alaska)\n\
rin(xy,time): [mm] double\n\
    Integrated rain over the last 24 hours.\n\
    (24h from solar noon to solar noon; or 2200UTC in Alaska)\n\
imonth(time): int (1-12)\n\
    Month of each timepoint.\n\
\n\
Optional Keyword Arguments\n\
--------------------------\n\
ffmc0: double\n\
    Initial value of fine fuel fuel moisture code (FMC) at the\n\
    beginning of the season.\n\
dmc0: double\n\
    Initial value of the duff moisture code (DMC) at the\n\
    beginning of the season.\n\
dc0: double\n\
    Initial value of the drought code (DC) at the\n\
    beginning of the season.\n\
debug: bool\n\
    Print additional debugging to stderr\n\
\n\
Returns: (bui, ffmc, isi, fwi, dmc, dc)\n\
---------------------------------------\n\
bui(xy,time): Buildup Index (double)\n\
    Numeric rating of the total\n\
    amount of fuel available for combustion. It is based on the\n\
    DMC and the DC. The BUI is generally less than twice the DMC\n\
    value, and moisture in the DMC layer is expected to help\n\
    prevent burning in material deeper down in the available fuel.\n\
\n\
ffmc(xy,time): Fine Fuel Moisture Code (double)\n\
    Numeric rating of the\n\
    moisture content of litter and other cured fine fuels. This\n\
    code is an indicator of the relative ease of ignition and the\n\
    flammability of fine fuel.\n\
\n\
isi(xy,time): Initial Spread Index (double)\n\
    Numeric rating of the expected rate of fire spread. It is\n\
    based on wind speed and FFMC. Like the rest of the FWI system\n\
    components, ISI does not take fuel type into account. Actual\n\
    spread rates vary between fuel types at the same ISI.\n\
\n\
fwi(xy,time): Fire Weather Index (double)\n\
    Numeric rating of fire intensity. It is based on the ISI and\n\
    the BUI, and is used as a general index of fire danger\n\
    throughout the forested areas of Canada.\n\
\n\
dmc(xy,time): Duff Moisture Code (double)\n\
    The Duff Moisture Code (DMC) is a numeric rating of the\n\
    average moisture content of loosely compacted organic layers\n\
    of moderate depth. This code gives an indication of fuel\n\
    consumption in moderate duff layers and medium-size woody\n\
    material.\n\
\n\
dc(xy,time): Drought Code (double)\n\
    The Drought Code (DC) is a numeric rating of the average\n\
    moisture content of deep, compact organic layers. This code is\n\
    a useful indicator of seasonal drought effects on forest fuels\n\
    and the amount of smoldering in deep duff layers and large\n\
    logs.";

static PyObject* canadafire_canadafire(PyObject *module, PyObject *args, PyObject *kwargs)
{
    // Input arrays
    PyArrayObject *tin;  // Temperature [C]
    PyArrayObject *hin;  // Relative Humidity [%, 0-100]
    PyArrayObject *win;  // Wind speed [km/h]
    PyArrayObject *rin;  // 24-hour integrated rain [mm]
    PyArrayObject *imonth = NULL;    // Month of each point in time
    int debug = 0;          // Optional kwarg, are we debugging?
    double ffmc00 = 85.5;    // Initial fine fuel moisture code (FMC)
    double dmc00 = 6.0;      // Initial duff moister code (DMC)
    double dc00 = 15.0;      // Initial drought code (DC)

    // List must include ALL arg names, including positional args
    static char *kwlist[] = {
        "tin", "hin", "win", "rin",    // *args
        "imonth", "ffmc0", "dmc0", "dc0", "debug",         // **kwargs
        NULL};
    // -------------------------------------------------------------------------    
    /* Parse the Python arguments into Numpy arrays */
    // p = "predicate" for bool: https://stackoverflow.com/questions/9316179/what-is-the-correct-way-to-pass-a-boolean-to-a-python-c-extension
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!|O!dddp",
        kwlist,
        &PyArray_Type, &tin,    // tin(yx,t)
        &PyArray_Type, &hin,    // hin(yx,t)
        &PyArray_Type, &win,    // win(yx,t)
        &PyArray_Type, &rin,    // rin(yx,t)
//        &PyArray_Type, &imonth,    // imonth(t)
        &PyArray_Type, &imonth, &ffmc00, &dmc00, &dc00, &debug
        )) return NULL;

    // Construct standard 183-day imonth array if none given.
    if (imonth == NULL) {
        static const int imonth_dims[] = {183};
        static const int imonth_strides[] = {sizeof(int)};
        static const int imonth_data[] = {
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,      // April
            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,    // May
            6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,      // June
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,    // July
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,    // August
            9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9};     // September
        imonth = (PyArrayObject*) PyArray_NewFromDescr(&PyArray_Type, 
            PyArray_DescrFromType(NPY_INT), 1, imonth_dims, imonth_strides, imonth_data,
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
    }




    const char *input_names[] = {"tin", "hin", "win", "rin"};
    PyArrayObject *inputs0[] = {tin, hin, win, rin};
    const int ninputs = 4;
    char msg[256];

    if (debug) {
        fprintf(stderr, "module = ");
        PyObject_Print(module, stderr, 0);
        fprintf(stderr, "\nargs = ");
        PyObject_Print(args, stderr, 0);
        fprintf(stderr, "\nkwargs = ");
        PyObject_Print(kwargs, stderr, 0);
        fprintf(stderr, "\n");
    }

    // -----------------------------------------
    // Make sure main arrays all have same rank and dimensions
    for (int i=0; i<ninputs; ++i) {
        if (!check_rank_dims(input_names[i], inputs0[i], input_names[0], inputs0[0]))
            return NULL;
    }

    // -----------------------------------------
    // Collapse space dimensions

    // Save original shape
    int ndim0 = PyArray_NDIM(tin);
    npy_intp const ntime = PyArray_DIM(tin, ndim0-1);
    npy_intp _dims0[ndim0];
    for (int j=0; j<ndim0; ++j) _dims0[j] = PyArray_DIM(tin,j);
    PyArray_Dims shape0 = {_dims0, ndim0};

    // Reshape to rank 2
    npy_intp nxy = 1;
    {
        for (int j=0; j<ndim0-1; ++j) nxy *= PyArray_DIM(inputs0[0],j);
        npy_intp _dims1[] = {nxy, ntime};
        PyArray_Dims shape1 = {_dims1, 2};

        // This will copy arrays if not already in C order.
        // That would be a good thing, since it would make time the lowest-stride dimension.
        tin = (PyArrayObject *)PyArray_Newshape(tin, &shape1, NPY_CORDER);
        if (!tin) return NULL;
        hin = (PyArrayObject *)PyArray_Newshape(hin, &shape1, NPY_CORDER);
        if (!hin) return NULL;
        win = (PyArrayObject *)PyArray_Newshape(win, &shape1, NPY_CORDER);
        if (!win) return NULL;
        rin = (PyArrayObject *)PyArray_Newshape(rin, &shape1, NPY_CORDER);
        if (!rin) return NULL;
    }

    // -------------------------------------------------
    // imonth: Check type, rank and size
    // Type
    if (PyArray_DESCR(imonth)->type_num != NPY_INT) {
        PyErr_SetString(PyExc_TypeError, "Parameter imonth must have type int");
        return NULL;
    }

    // Rank
    if (PyArray_NDIM(imonth) != 1) {
        PyErr_SetString(PyExc_TypeError,
            "Parameter imonth must have rank 1");
        return NULL;
    }

    // Dimensions
    if ((int)PyArray_DIM(imonth, 0) != (int)ntime) {
        sprintf(msg, "Parameter imonth must have length %d equal to the time dimension of other variables; but has length %d instead.", (int)ntime, (int)PyArray_DIM(imonth,0));
        PyErr_SetString(PyExc_TypeError, msg);
        return NULL;
    }

    // -------------------------------------------------------------------------    

    // Output arrays
    PyArrayObject *ffmcout = (PyArrayObject *)PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);
    PyArrayObject *isiout = (PyArrayObject *) PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);
    PyArrayObject *fwiout = (PyArrayObject *)PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);
    PyArrayObject *dmcout = (PyArrayObject *)PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);
    PyArrayObject *buiout = (PyArrayObject *)PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);
    PyArrayObject *dcout = (PyArrayObject *)PyArray_NewLikeArray(tin, NPY_ANYORDER, NULL, 0);


    // Loop through space; typically a 2D grid.  But could be (for
    // example) a small set of discrete regions for which we have
    // appropriate input data.
    for (int ii=0; ii<nxy; ++ii) {

        // enter in the inital values here
        // 3 fuel moisture codes, initialize
        double ffmc0 = ffmc00;
        double dmc0 = dmc00;
        double dc0 = dc00;


        // Loop through time.  Assumed to be solar noon on each of
        // April 1 -- September 30.  But it could be any set of
        // timepoints.
        for (int ti=0; ti<ntime; ++ti) {

            // Obtain values at this point in (space,time).
            const double t = *((double *)PyArray_GETPTR2(tin,ii,ti));
            const double h = *((double *)PyArray_GETPTR2(hin,ii,ti));
            const double w = *((double *)PyArray_GETPTR2(win,ii,ti));
            const double r = *((double *)PyArray_GETPTR2(rin,ii,ti));
            const int im =  *((int *)PyArray_GETPTR1(imonth,ti));

            // Run the core computation on a single gridpoint.
            double bui, ffm, isi, fwi, dsr, dmc, dc;   // Output vars
            canadafire(
                t, h, w, r, im,
                ffmc0, dmc0, dc0,    // Running values from one day to the next
                &bui, &ffm, &isi, &fwi, &dsr, &dmc, &dc);

            // housekeeping items
            // set todays values to the yesterdays values before going on
            ffmc0 = ffm;
            dmc0 = dmc;
            dc0 = dc;

            // Store results back in Arrays
            *((double *)PyArray_GETPTR2(buiout,ii,ti)) = bui;
            *((double *)PyArray_GETPTR2(ffmcout,ii,ti)) = ffm;
            *((double *)PyArray_GETPTR2(isiout,ii,ti)) = isi;
            *((double *)PyArray_GETPTR2(fwiout,ii,ti)) = fwi;
            *((double *)PyArray_GETPTR2(dmcout,ii,ti)) = dmc;
            *((double *)PyArray_GETPTR2(dcout,ii,ti)) = dc;
            
        }
    }

    // ---------------------------------------------------------
    // Reshape to original rank
    ffmcout = (PyArrayObject *)PyArray_Newshape(ffmcout, &shape0, NPY_ANYORDER);
    if (!ffmcout) return NULL;
    isiout = (PyArrayObject *)PyArray_Newshape(isiout, &shape0, NPY_ANYORDER);
    if (!isiout) return NULL;
    fwiout = (PyArrayObject *)PyArray_Newshape(fwiout, &shape0, NPY_ANYORDER);
    if (!fwiout) return NULL;
    dmcout = (PyArrayObject *)PyArray_Newshape(dmcout, &shape0, NPY_ANYORDER);
    if (!dmcout) return NULL;
    buiout = (PyArrayObject *)PyArray_Newshape(buiout, &shape0, NPY_ANYORDER);
    if (!buiout) return NULL;
    dcout = (PyArrayObject *)PyArray_Newshape(dcout, &shape0, NPY_ANYORDER);
    if (!dcout) return NULL;

    // Return a tuple of the output arrays we created.
    // https://stackoverflow.com/questions/3498210/returning-a-tuple-of-multipe-objects-in-python-c-api
    return PyTuple_Pack(6, buiout,ffmcout,isiout,fwiout,dmcout,dcout);
}
// ============================================================
// Random other Python C Extension Stuff
static PyMethodDef CanadafireMethods[] = {
    {"canadafire",
        canadafire_canadafire,
        METH_VARARGS | METH_KEYWORDS, canadafire_canadafire_docstring},
    {NULL, NULL, 0, NULL}
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "canadafire",    // Name of module
    module_docstring,    // Per-module docstring
    -1,  /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    CanadafireMethods,    // Functions
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_canadafire(void)
{
    import_array();    // Needed for Numpy

    return PyModule_Create(&moduledef);
}
