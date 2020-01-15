# User Defined Function for Droplet Phase Detection for VOF Model in ANSYS Fluent (cpad_udf_library)
UDF library for the detection and localization of separated continuos phase areas (for example droplets) in Multiphase VoF (Volume of Fluid) ANSYS Fluent cases.

Currently Limitations:
  - limited to serial double precision Fluent cases

Settings:
``` cpp
#define _PHASE_IDX 0            /* for phase to detect droplets of */
#define _MIN_VOL_FRAC 0.01      /* Lower Limit for phase detection */
#define _DETECT_OVER_EDGES 1    /* 0 = detect over faces */

#define _FLUID_  1              /* Fluid Cell Zone ID*/
```


Captured properties of separated continuos phase areas are:
  - cell count
  - x,y,z center of mass
  - volume
  - mass
  - average volume fraction
  - max volume fraction
  - connected boundary id's