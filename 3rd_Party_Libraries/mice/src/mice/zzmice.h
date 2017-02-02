/*

-Abstract

   The structure and protype definitions corresponding to the argument
   check code in zzmice.c.

   Do not edit, touch, or in any way alter the code within this file.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   MICE.REQ

-Keywords

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.
   
-Parameters
*/

#ifndef zzmice_h
#define zzmice_h


/*
-Exceptions

   None.

-Files

   None.

-Particulars

   None.

-Examples

   None.

-Restrictions

   None.

-Literature_References

   MICE.REQ

-Author_and_Institution

   E.D. Wright        (JPL)
   G. Chinn           (JPL)

-Version

   -Mice Version 1.0.2 31-OCT-2012 (EDW)
   
      Added MiceFrinfo, MiceSFS, MicePVN enums.

   -Mice Version 1.0.1 17-DEC-2008 (EDW)

      Addition of MiceWnsumd enum.

   -Mice Version 1.0.0 14-FEB-2008 (EDW)

-Index_Entries

*/


/* 

The structure defining the argument check information, ArgCheck.

-Examples

   Declaration of ArgCheck for the subsol_c interface.
   
   struct argcheck ArgCheck[] = 
      {
      { "method", MiceChar  , 0, {0}, 0},
      { "target", MiceChar  , 0, {0}, 0},
      { "et"    , MiceDouble, 0, {0}, 1},
      { "abcorr", MiceChar  , 0, {0}, 0},
      { "obsrvr", MiceChar  , 0, {0}, 0},
      { "spoint", MiceDouble, 1, {3}, 1},
      };

*/

/*
Enums used to tag argument data types in the ArgCheck arrays. Use 
"MiceIgnore"  to prevent allocation to an output pointer (plhs), when 
used, the interface needs to peform memory allocation to the pointer.
*/
enum MiceType
   {
   MiceInvalidType,
   MiceNameID,
   MiceFrinfo,
   MiceState,
   MicePos,
   MiceNear,
   MiceSurf,
   MicePlane,
   MiceEllipse,
   MicePool,
   MiceChar,
   MiceDouble,
   MiceInt,
   MiceBoolean,
   MiceWin,
   MiceSub_PS,
   MiceSurf_PS,
   MiceIlum,
   MiceWnsumd,
   MiceSFS,
   MicePVN,
   MiceIgnore,
   };


struct argcheck {
                char             * name;      /* Variable name. */
                          
                enum MiceType      type;      /* 
                                                 Variable type identifier,
                                                 an enum. 
                                              */
                                              
                int            min_dims;      /* 
                                                 The dimensinality of the 
                                                 variable. 
   
                                                    0 for a scalar
                                                    1 for a Nx1 array, 
                                                    2 for a NxM array.
                                              */
                                                        
                SpiceInt           dims[4];   /* The expected dimension and 
                                                 corresponding size of each
                                                 dimension. 
                                              */
                                                        
                SpiceInt   is_vectorizable;   /* Flag to mark if the variable
                                                 may pass in a vectorized state.
                                    
                                                    1 to indicate yes
                                                    0 to indicate no
                                              */
                };


/* 

The structure defining the vectorization state. Returned by 'mice_checkargs'.
 
*/
struct           extra_dims {
                            SpiceInt  count;
                            SpiceInt  first_vector_arg_index;
                            SpiceInt  vectorized [32];
                            SpiceInt  offset     [32];
                            };


void           check_arg_num( int x_nrhs, int x_nlhs, int nrhs, int nlhs );

void           mice_fail( long cnt );

void           struct_fields( enum   MiceType        type , 
                             SpiceInt              * n    , 
                             const  char         *** names,
                             const  enum MiceType ** types,
                             const       int      ** sizes );

struct         extra_dims  * mice_checkargs(int                 nlhs,
                                            mxArray           * plhs[],
                                            int                 nrhs,
                                            const mxArray     * prhs[],
                                            struct argcheck   * argcheck);


#endif

