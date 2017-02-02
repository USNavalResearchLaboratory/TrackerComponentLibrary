/*
-Procedure zzmice_CreateIntScalar ( Scalar integer to MATLAB scalar integer)

-Abstract

   Simple function to copy an integer (SpiceInt) to a MATLAB scalar integer,
   i.e. a 1x1 integer array.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
   KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   MICE.REQ

-Keywords

   MATLAB

-Brief_I/O

   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   n          I   Integer value for copy to a MATLAB scalar
   mx_n       O   Pointer to an mxArray containing 'n'

-Detailed_Input

   n          a scalar integer to copy to a MATLAB form of scalar integer

-Detailed_Output

   mx_n       a pointer to the new mxArray, defined as a 1x1 array, 
              containing 'n'

-Parameters

   None.

-Exceptions

   None.

-Files

   None.

-Particulars

   This function exists only to wrap the two calls required to create a
   MATLAB scalar integer.

   Mice functionality depends on a number of routines provided with MATLAB by
   Mathworks. These routines (descriptions quoted from mathworks.com):

   -mxCreateNumericArray returns:
   
      A pointer to the created mxArray, if successful. If unsuccessful in a
      stand-alone (non MEX-file) application, mxCreateNumericArray returns NULL.
      If unsuccessful in a MEX-file, the MEX-file terminates and control returns
      to the MATLAB prompt. mxCreateNumericArray is unsuccessful when there is
      not enough free heap space to create the mxArray.

-Examples

   /.
   Assign 'n' to struct field "dim".
   ./
   mxSetField( plhs[0], i,"dim", zzmice_CreateIntScalar(n) );

-Restrictions

   See particulars.

-Literature_References

   Mathworks External Interfaces Reference:
   
      http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/

-Author_and_Institution

   E. D. Wright    (JPL)

-Version

   -Mice Version 1.0.0, 30-MAR-2007 (EDW)

-Index_Entries

   MATLAB

-&
*/

#include <stdlib.h>
#include "mex.h"
#include "SpiceUsr.h"

mxArray * zzmice_CreateIntScalar( SpiceInt n)
   {

   int                     sizearray[] = {1,1};
   mxArray               * mx_n;

   /*
   Rather a lot of code to copy a scalar.
   */
   mx_n = mxCreateNumericArray( 2, sizearray, mxINT32_CLASS, mxREAL);

   *(SpiceInt*)mxGetData(mx_n) = n;

   return mx_n;
   }


