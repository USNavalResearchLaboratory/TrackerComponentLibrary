gpu = false; % compile the gpu version of SCS
float = false; % using single precision (rather than double) floating points
int = false; % use 32 bit integers for indexing
% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
% MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t
flags.LCFLAG = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX';
flags.INCS = '';
flags.LOCS = '';

%The path has been modified by David F. Crouse to reflect the folder
%structure for inclusion in the Tracker Component Library, January 2018.
common_scs = '../src/linalg.c ../src/cones.c ../src/accel.c ../src/cs.c ../src/util.c ../src/scs.c ../src/ctrlc.c ../src/normalize.c ../src/scs_version.c ../linsys/common.c scs_mex.c';
if (~isempty (strfind (computer, '64')))
    flags.arr = '-largeArrayDims';
else
    flags.arr = '';
end

if ( isunix && ~ismac )
    flags.link = '-lm -lut -lrt';
elseif  ( ismac )
    flags.link = '-lm -lut';
else
    flags.link = '-lut';
    flags.LCFLAG = sprintf('-DNOBLASSUFFIX %s', flags.LCFLAG);
end

if (float)
    flags.LCFLAG = sprintf('-DSFLOAT %s', flags.LCFLAG);
end
if (int)
    flags.INT = '';
else
    flags.INT = '-DDLONG';
end

if (flags.COMPILE_WITH_OPENMP)
    flags.link = strcat(flags.link, ' -lgomp');
end

%compile_direct and the GPU option commented out by David F. Crouse for
%inclusion in the Tracker Component Library, January 2018.
%compile_direct(flags, common_scs);
compile_indirect(flags, common_scs);
%if (gpu)
%    compile_gpu(flags, common_scs);
%end

%The path has been modified by David F. Crouse to reflect the folder
%structure for inclusion in the Tracker Component Library, January 2018.
% compile scs_version
mex -O -I../include ../src/scs_version.c scs_version_mex.c -output scs_version 

%Test run commented out by David F. Crouse for inclusion in the Tracker
%Component Library, January 2018. Some older versions of Matlab and do not
%register the new mex files until the script ends, meaning that the test
%run would always fail.

% clear data cones
% disp('Example run:');
% randn('seed',9)
% m = 9;
% n = 3;
% data.A = sparse(randn(m,n));
% data.b = randn(m,1);
% data.c = randn(n,1);
% cones.l = m;
% [x,y,s,info] = scs_indirect(data,cones,[]);
% [x,y,s,info] = scs_direct(data,cones,[]);
% 
% 
% if (gpu)
%     [x,y,s,info] = scs_gpu(data,cones,[]);
% end
% 
% % test-warm start with solution
% disp('Warm-starting:')
% data.x = x;
% data.y = y;
% data.s = s;
% [x,y,s,info] = scs_indirect(data,cones,[]);
% 
% 
% disp('SUCCESSFULLY INSTALLED SCS')
% disp('(If using SCS with CVX, note that SCS only supports CVX v3.0 or later).')
