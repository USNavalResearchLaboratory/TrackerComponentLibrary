This code is taken from version 2.6.1 of the NLOpt library from
http://ab-initio.mit.edu/wiki/index.php/NLopt
The code has been modified to get rid of most of the warning when 
compiling. These include, for example, complaints that variables are unused
or are used prior to being initialized. Also the memory allocation and
freeing routines in the library have been changed to use the methods from
Matlab as opposed to those in stdlib. Additionally, the external
dependencies that allow the user to limit the execution time have been
removed. An additional input has also been added to the direct_optimize
function to allow a user to change the maximum number of rectangle
divisions. Previously, it had been hard-coded in the function. Changes that
are not meant to silence warnings are commented with DFC in the relevant
files. The DIRparallel file is not used in the Tracker Component Library.

August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
