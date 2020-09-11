# scs-matlab
Matlab interface for [SCS](https://github.com/cvxgrp/scs) 2.0.0 and higher.

To install SCS in Matlab from source:
```shell
git clone --recursive https://github.com/bodono/scs-matlab.git
```
Then in a Matlab session
```matlab
cd <path/to/scs-matlab>
make_scs
```
This will produce two mex files, one for the direct version of SCS and one for
the indirect version.

Remember to include the `scs-matlab` directory in your Matlab path if you wish
to use the mex file in your Matlab code. The calling sequence is (for the
indirect version):
```matlab
[x,y,s,info] = scs_indirect(data, cones, settings)
```

where data is a struct containing `A`, `b`, and `c`, settings is a struct
containing solver options (see matlab file, can be empty),
and cones is a struct that contains one or more of:
+ `f`  (num primal zero / dual free cones, i.e. primal equality constraints)
+ `l`  (num linear cones)
+ `q`  (array of SOCs sizes)
+ `s`  (array of SDCs sizes)
+ `ep` (num primal exponential cones)
+ `ed` (num dual exponential cones)
+ `p`  (array of primal/dual power params).

Type `help scs_indirect` at the Matlab prompt to see its documentation.

