/*

-Procedure zzmice ( Mice auxiliary routines )

-Abstract

   None.

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

   None.

-Keywords

   None.

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters

   None.

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

   None.

-Author_and_Institution

   E. D. Wright    (JPL)

-Version

   -Mice Version 1.2.0 30-OCT-2012 (EDW)

      Added definitions and conditions for MiceFrinfo, MiceSFS,
      MicePVN structures.

      Structure MiceBodID renamed MiceNameID.

      Structure MiceSurf renamed MiceNear.

   -Mice Version 1.1.0 15-DEC-2008 (EDW)

      Added definitions and conditions for MicePlane, MiceEllipse,
      and MiceWnsumd structures.

   -Mice Version 1.0.0 12-FEB-2008 (EDW)

-Index_Entries

   MATLAB

-&
*/

#include <stdio.h>
#include "SpiceUsr.h"
#include "SpiceZmc.h"
#include "mex.h"
#include "mice.h"
#include "zzmice.h"
#include "zzerror.h"


/*

   Define structures used as Mice return types. Current code limited to
   twelve return fields.

         MiceNameID:
            bodc2n_c
            bodn2c_c
            bods2c_c
            cnmfrm_c

         MiceFrinfo:
            frinfo_c

         MiceState:
            spkezr_c

         MicePos:
            spkpos_c

         MiceNear:
            subpt_c
            nearpt_c
            npedln_c

         MiceSurf:
            srfxpt_c

         MicePlane:
            inedpl_c
            inelpl_c
            inrypl_c
            nvc2pl_c
            nvp2pl_c
            pjelpl_c
            pl2nvc_c
            pl2nvp_c
            pl2psv_c
            psv2pl_c
            vprjp_c
            vprjpi_c

         MiceEllipse:
            cgv2el_c
            edlimb_c
            el2cgv_c
            inedpl_c
            inelpl_c
            npelpt_c
            pjelpl_c

         MicePool:
            dtpool_c

         MiceSub_PS:
            subpnt_c
            subslr_c

         MiceSurf_PS:
            sincpt_c

         MiceIlum:
            ilumin_c

         MiceWnsumd:
            wnsumd_c

         MiceSFS:
            spksfs_c

         MicePVN:
            spkpvn_c

*/
static struct structdata
   {
   enum  MiceType    type;
   const char      * name;
   int               nFields;
   const char      * fields    [12];
   enum  MiceType    fieldTypes[12];
   int               fieldSizes[12];
   }
   structdata[] =
      {
         {
           MiceNameID,
           "MiceNameID",
           3,
           {"code",  "name",   "found"    },
           {MiceInt, MiceChar, MiceBoolean},
           {1,       1,        1          }
         },

         {
           MiceFrinfo,
           "MiceFrinfo",
           4,
           {"center", "class", "class_ID", "found"    },
           {MiceInt,  MiceInt, MiceInt,    MiceBoolean},
           {1,        1,        1,         1          }
         },

         {
           MiceState,
           "MiceState",
           2,
           {"state",    "lt"      },
           {MiceDouble, MiceDouble},
           {6,          1         }
         },

         {
           MicePos,
           "MicePos",
           2,
           {"pos",      "lt"      },
           {MiceDouble, MiceDouble},
           {3,          1         }
         },

         {
           MiceNear,
           "MiceNear",
           2,
           {"pos",      "alt"     },
           {MiceDouble, MiceDouble},
           {3,          1         }
         },

         {
           MiceSurf,
           "MiceSurf",
           5,
           {"spoint",   "dist",     "trgepc",   "obspos",   "found"    },
           {MiceDouble, MiceDouble, MiceDouble, MiceDouble, MiceBoolean},
           {3,          1,          1,          3,          1          }
         },

         {
           MicePlane ,
           "MicePlane",
           2,
           {"normal",   "constant"},
           {MiceDouble, MiceDouble},
           {3,          1         }
         },

         {
           MiceEllipse,
           "MiceEllipse",
           3,
           {"center",   "semiMajor", "semiMinor"},
           {MiceDouble, MiceDouble,  MiceDouble },
           {3,          3,           3          }
         },

         {
           MicePool  ,
           "MicePool",
           3,
           {"n",        "type",   "found"    },
           {MiceDouble, MiceChar, MiceBoolean},
           {1,          1,        1          }
         },

         {
           MiceSub_PS ,
           "MiceSub_PS",
           3,
           {"spoint",   "trgepc",    "srfvec"  },
           {MiceDouble, MiceDouble,  MiceDouble},
           {3,          1,           3         }
         },

         {
           MiceSurf_PS ,
           "MiceSurf_PS",
           4,
           {"spoint",   "trgepc",   "srfvec",   "found"    },
           {MiceDouble, MiceDouble, MiceDouble, MiceBoolean},
           {3,          1,          3,          1          }
         },

         {
           MiceIlum ,
           "MiceIlum",
           5,
           {"trgepc",   "srfvec",   "phase",   "solar",     "emissn"  },
           {MiceDouble, MiceDouble, MiceDouble, MiceDouble, MiceDouble},
           {1,          3,          1,          1,          1         }
         },

         {
           MiceWnsumd,
           "MiceWnsumd",
           5,
           {"meas",     "avg",      "stddev",   "shortest", "longest"},
           {MiceDouble, MiceDouble, MiceDouble, MiceInt,    MiceInt  },
           {1,          1,          1,          1,          1        }
         },

         {
           MiceSFS ,
           "MiceSFS",
           4,
           {"handle", "descr",   "ident",   "found"    },
           {MiceInt,  MiceDouble, MiceChar, MiceBoolean},
           {1,        5,          1,        1          }
         },

         {
           MicePVN,
           "MicePVN",
           3,
           {"ref",    "state",    "center"},
           {MiceInt,  MiceDouble, MiceInt },
           {1,        6,          1       }
         },

      };

static int nStructdata = sizeof(structdata)/sizeof(structdata[0]);


/*

-Procedure check_arg_num  ( Check number of arguments in a Mice interface call )

-Abstract

   Test the number of expected input and output arguments match the
   number of arguments sent from MATLAB. Error out if the values
   do not conform to the expected values.

   This error indicates an error in the code of the .m file used to call
   a Mice interface.

-Disclaimer

   None.

-Required_Reading

   None.

-Keywords

   None.

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters

   None.

-Exceptions

   None.

-Files

   None.

-Particulars

   The .m wrapper (should) check for improper argument lists from the user.
   This function checks for improper argument lists implemented by the
   system designer.

-Examples

   None.

-Restrictions

   None.

-Literature_References

   None.

-Author_and_Institution

   None.

-Version

   None.

-Index_Entries

   None.

-&

*/
void check_arg_num ( int x_nrhs, int x_nlhs, int nrhs, int nlhs )
   {

   /*
   A short string for the error message.
   */
   SpiceChar   message[128];

   /*
   The expected right-hand argument list count includes the called function
   name at index 0. Add one to the Mice defined number of arguments to
   compensate.
   */
   if ( x_nrhs != (nrhs+1) )
      {
      sprintf( message,
         "This routine requires %d input argument(s). %d inputs assigned",
         nrhs,
         (x_nrhs-1) );

      mexErrMsgTxt( message );
      }
   else if ( x_nlhs != nlhs )
      {
      sprintf( message,
         "This routine requires %d output argument(s). %d outputs assigned",
          nlhs,
          x_nlhs );

      mexErrMsgTxt( message );
      }

   }



/*

-Procedure mice_fail (report SPICE errors)

-Abstract

   Respond to SPICE errors by building the traceback string, passing that string
   and the error message to MATLAB, then reset the SPICE error system.
   the zzerror() call performs all error subsystem calls.

-Disclaimer

   None.

-Required_Reading

   None.

-Keywords

   None.

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters

   None.

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

   None.

-Author_and_Institution

   None.

-Version

   None.

-Index_Entries

   None.

-&

*/
void mice_fail( long cnt )
   {

   /*
   Count inputs from the interface routine, a C routine so keyed to
   the convention of zero based arrays. Add INDEX_BASE to the count
   to compensate since the index count should reflect a one base count.
   */
   mexErrMsgTxt( zzerror( cnt + INDEX_BASE) );
   }



/*
 returns an enum flag corresponding to a particular cspice library
 structure type given an mxArray which contains a struct.
 this routine assumes that the mxArray actually contains a struct!!!
 (not a double/char)
*/
void struct_fields( enum   MiceType          type,
                           SpiceInt        * n    ,
                    const  char          *** names,
                    const  enum MiceType  ** types,
                    const  int            ** sizes)
   {
   int                     i;

   for (i=0;i<nStructdata;++i)
      {
      if (type == structdata[i].type)
         {
         *n     = structdata[i].nFields;
         *names = (const char**        )structdata[i].fields;
         *types = (const enum MiceType*)structdata[i].fieldTypes;
         *sizes = (const int*          )structdata[i].fieldSizes;
         return;
         }
      }

   return;
   }



/*
-Procedure mice_checkargs ( argument checking routine )

-Abstract

   None.

-Disclaimer

   None.

-Required_Reading

   None.

-Keywords

   None.

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters

   None.

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

   None.

-Author_and_Institution

   None.

-Version

   None.

-Index_Entries

   None.

-&
*/
struct extra_dims  * mice_checkargs(int                 nlhs,
                                    mxArray           * plhs[],
                                    int                 nrhs,
                                    const mxArray     * prhs[],
                                    struct argcheck   * argcheck)
   {
   SpiceInt                    ii;
   SpiceInt                i;
   SpiceInt                j;
   SpiceInt                k;
   SpiceInt                m          = 0;
   SpiceInt                n          = 0;
   SpiceInt                arg_dimension;
   SpiceInt                num_elem;
   SpiceInt                nfields;
   SpiceInt                vec_flag    = 0;
   SpiceInt                vec_measure = 0;
   SpiceInt                is_a_struct;
   int                     sizearray[3];
   const int             * Dims;
   const int             * fsizes;

   char                    msg[1024];
   const char           ** fnames;

   const enum MiceType   * ftypes;

   /* for vectorizable functions */
   static struct           extra_dims extra;

   /*
   Initialize  the 'extra' fields.
   */
   extra.count = 0;

   /*
   Clunky implementation, but obvious as regards function.
   */
   for (i=0; i<32; i++)
      {
      extra.vectorized[i] = -1;
      extra.offset[i]     =  0;
      }


   /*
   Argument check index for descriptions in ArgCheck.
   */
   ii = 0;

   /*

   INPUT - Loop over the list of input arguments.

   vec_flag is a flag which will indicate whether the output is
   vectorized.

   if the function is vectorizable then:

     vec_flag = 1 if any of the input variables is a vector.
     vec_flag = 0 if all of the input variables are scalars.

   If the function is not vectorizable then vec_flag = 0.
   */

   for (i=1; i<nrhs; i++)
      {

      /*
      Confirm the input has a known type.
      */
      if (mxGetClassID(prhs[i]) != mxDOUBLE_CLASS  &&
          mxGetClassID(prhs[i]) != mxINT32_CLASS   &&
          mxGetClassID(prhs[i]) != mxCHAR_CLASS    &&
          mxGetClassID(prhs[i]) != mxSTRUCT_CLASS
         )
         {
         sprintf( msg,
                  "MICE(BUG): Input argument (`%s') is "
                  "not a recognized Class type. "
                  "This indicates an interface bug. Contact NAIF.",
                  argcheck[ii].name);

         mexErrMsgTxt(msg);
         }



      /*
      Check the type of the input argument against the expected type.
      Signal an error for cases of type mismatch.
      */

      switch ( argcheck[ii].type )
         {

         case MiceDouble:
         case MiceWin:

            /*
            All elements of the argument should have type double precision.
            */

            if ( mxGetClassID(prhs[i]) != mxDOUBLE_CLASS )
               {
               sprintf(msg,
                  "MICE(BADVAL): All elements of input argument (`%s') "
                  "must have type double.",
                  argcheck[ii].name);

               mexErrMsgTxt(msg);
               }

            break;

         case MiceInt:

            /*
            Argument should have type integer.
            */

            if (  mxGetClassID(prhs[i]) != mxINT32_CLASS )
               {
               sprintf(msg,
                  "MICE(BADARG): Input argument (`%s') "
                  "must be an integer value.",
                  argcheck[ii].name);

               mexErrMsgTxt(msg);
               }

            break;

         case MiceBoolean:

            /*
            Argument should have type integer... the interface wrapper
            converts MATLAB booleans to integer. This should never signal.
            */

            if (  mxGetClassID(prhs[i]) != mxINT32_CLASS )
               {
               sprintf(msg,
                  "MICE(BADARG): Input argument (`%s') "
                  "must be an integer representation of a boolean value."
                  "Check the interface wrapper file for a logical to "
                  "int conversion.",
                  argcheck[ii].name);

               mexErrMsgTxt(msg);
               }

            break;

         case MiceChar:

            /*
            Argument should have type character (string or single char).
            */

            if (mxGetClassID(prhs[i]) != mxCHAR_CLASS)
               {
               sprintf(msg,
                  "MICE(BADARG): Input argument (`%s') must be alphabetic.",
                  argcheck[ii].name);
               mexErrMsgTxt(msg);
               }

            break;

         case MiceEllipse:
         case MicePlane:

            /*
            All elements of the argument should have type structure.
            */

            if ( mxGetClassID(prhs[i]) != mxSTRUCT_CLASS )
               {
               sprintf(msg,
                  "MICE(BADVAL): All elements of input argument (`%s') "
                  "must have type struct.",
                  argcheck[ii].name);

               mexErrMsgTxt(msg);
               }

            break;

         default:

            /*
            Unknown argument type. This error indicates a mistake when coding a
            new argument type.
            */
            sprintf( msg,
                     "MICE(BUG): Input argument (`%s') of unknown type."
                     "This occurs only if mice.h does not list a argument "
                     "type coded into one of the interface's ArgCheck "
                     "assignment.",
                     argcheck[ii].name);

            mexErrMsgTxt(msg);

            break;
         }


      /*
      Determine input argument type, scalar or vector, from the size of
      each input argument's dimension. All non vectorized matrices input
      arguments must have dimension 2, regardless of whether the argument
      contains a scalar [1x1], a vectorized scalar N-row-vector [1xN], an
      N-column-vector [Nx1], or a string array [Nxlength_of_string]
      (a special case).

      Retrieve the dimensionality (number of degrees of freedom) of
      the argument.
      */
      arg_dimension = mxGetNumberOfDimensions(prhs[i]);

      /*
      Any argument to an interface must have either 2 or 3 dimensions.
      */
      if ( arg_dimension != 2 && arg_dimension != 3 )
         {
         sprintf(msg,
                "MICE(BUG): Incorrect number of dimensions for input argument "
                "(`%s'). Correct value either 2 or 3. Number of dimensions"
                "of this argument: %ld. This indicates an interface bug. "
                "Contact NAIF.",
                argcheck[ii].name,
                (long) arg_dimension);

         mexErrMsgTxt(msg);
         }


      if ( argcheck[ii].min_dims != 0 &&
           argcheck[ii].min_dims != 1 &&
           argcheck[ii].min_dims != 2 )
         {
         sprintf( msg,
                  "MICE(BUG): The min_dims value %d in argument `%s' not 0, "
                  "1, or 2. Check the argument specification in ArgCheck. "
                  "This indicates an interface bug. Contact NAIF.",
                  argcheck[ii].min_dims,
                  argcheck[ii].name);

         mexErrMsgTxt(msg);
         }



      /*
      Retrieve the dimension size array. Dims has length arg_dimension,
      each element of Dim containing the number of elements corresponding
      to that degree of freedom.
      */
      Dims     = mxGetDimensions(prhs[i]);
      num_elem = mxGetNumberOfElements(prhs[i]);

      /*
      Does this variable argument consist of double precision or integer value?
      */
      if ( mxGetClassID(prhs[i]) == mxDOUBLE_CLASS ||
           mxGetClassID(prhs[i]) == mxINT32_CLASS )
         {

         extra.vectorized[ii] = 0;

         /*
         The following blocks assumes error checks confirmed 'arg_dimension'
         has value 2 or 3.
         */
         if ( arg_dimension == 3 )
            {
            extra.vectorized[ii] = 1;
            }
         else
            {
            if(    Dims[0] == 1
                && Dims[1]  > 1
                && argcheck[ii].min_dims == 0 )
               {
               extra.vectorized[ii] = 1;
               }

            if(     (Dims[0]               == argcheck[ii].dims[0])
                &&  (argcheck[ii].min_dims == 1)
                &&  (Dims[1]               >  1) )
               {
               extra.vectorized[ii] = 1;
               }

            }

         /*
         If the argument might be vectorized, calculate the vectorization
         measure.
         */
         if ( argcheck[ii].is_vectorizable )
            {

            vec_measure = 0;

            switch ( argcheck[ii].min_dims )
               {
               case 0:
               case 1:
                  vec_measure = Dims[arg_dimension - 1];
                  break;

               case 2:
                  if( arg_dimension == 2 )
                     {
                     vec_measure = 1;
                     }
                  else
                     {
                     vec_measure = Dims[arg_dimension - 1];
                     }

                  break;

               default:
                  sprintf( msg,
                     "MICE(BUG): The min_dims value %d in argument `%s' not 0, "
                     "1, or 2. Check the argument specification in ArgCheck. "
                     "This indicates an interface bug. Contact NAIF.",
                     argcheck[ii].min_dims,
                     argcheck[ii].name);

                  mexErrMsgTxt(msg);
                  break;
               }

            /*
            Store the required measure of vectorization for all vectorized
            arguments.
            */
            if (extra.count == 0)
               {
               extra.first_vector_arg_index = ii;
               extra.count                  = vec_measure;
               }
            else if (extra.count != vec_measure)
               {
               sprintf( msg,
                        "MICE(BADARG): Input argument (`%s') "
                        "has vectorization measure %ld. It must "
                        "have same measure as `%s', %ld",
                        argcheck[ii].name,
                        (long) vec_measure,
                        argcheck[extra.first_vector_arg_index].name,
                        (long) extra.count
                        );

               mexErrMsgTxt(msg);
               }

            }


         if (!argcheck[ii].is_vectorizable && extra.vectorized[ii])
            {
            sprintf(msg,
               "Input argument (`%s') is not vectorizable.",
               argcheck[ii].name);

            mexErrMsgTxt(msg);
            }


         if ( argcheck[ii].is_vectorizable && extra.vectorized[ii] )
            {
            n = 2;

            switch ( argcheck[ii].min_dims )
               {
               case 0:
                  extra.offset[ii] = 1;
                  break;

               case 1:
                  extra.offset[ii] = argcheck[ii].dims[0];
                  break;

               case 2:
                  extra.offset[ii] = argcheck[ii].dims[0]
                                   * argcheck[ii].dims[1];
                  break;

               default:
                  sprintf( msg,
                     "MICE(BUG): The min_dims value %d in argument `%s' not 0, "
                     "1, or 2. Check the argument specification in ArgCheck. "
                     "This indicates an interface bug. Contact NAIF.",
                      argcheck[ii].min_dims,
                      argcheck[ii].name);

                  mexErrMsgTxt(msg);
                  break;
               }

            vec_flag = 1;

            }
         else
            {

            switch ( argcheck[ii].min_dims )
               {

               case 0:

                  if ( Dims[0]!=1 || Dims[1]!=1 )
                     {
                     sprintf(msg,
                        "Input argument (`%s') must be a "
                        "scalar (1x1) or a vectorized scalar (1xN).",
                        argcheck[ii].name);

                     mexErrMsgTxt(msg);
                     }

                  break;

               case 1:

                  m = argcheck[ii].dims[0];

                  if ( m != 0 )
                     {

                     /*
                     The vector length explicitly stated in ArgCheck.
                     */
                     if ( Dims[0]!=m || Dims[1]!=1 )
                        {
                        sprintf( msg,
                           "MICE(BADARG): Input argument (`%s') must be "
                           "an %ldx1 vector or a vectorized vector (%ldxN).",
                           argcheck[ii].name,
                           (long) m,
                           (long) m);

                        mexErrMsgTxt(msg);
                        }

                     }
                  else
                     {

                     /*
                     Confirm the input as an Nx1 array.
                     */
                     if (Dims[1]!=1)
                        {
                        sprintf( msg,
                           "Input argument (`%s') must be an "
                           "Nx1 array (column vector).",
                           argcheck[ii].name );

                        mexErrMsgTxt(msg);
                        }

                     }

                  break;

               case 2:

                  m = argcheck[ii].dims[0];
                  n = argcheck[ii].dims[1];

                  if (m !=0 && n !=0 )
                     {

                     /*
                     The ArgCheck block explicitly states array dimensions.
                     The input argument should have those dimensions.
                     */
                     if (Dims[0]!=m || Dims[1]!=n)
                        {
                        sprintf( msg,
                           "MICE(BADARG): Input argument (`%s') must "
                           "be an %ldx%ld matrix.",
                           argcheck[ii].name,
                           (long) m,
                           (long) n);

                        mexErrMsgTxt(msg);
                        }

                     }
                  else if ( m != 0 )
                     {

                     /*
                     N not stated in declaration:

                     min_dims = 2
                     dims = { M, 0}

                     */
                     if (Dims[0]!=m)
                        {
                        sprintf( msg,
                           "Input argument (`%s') must be a %ldXN matrix.",
                           argcheck[ii].name,
                           (long) m);

                        mexErrMsgTxt(msg);
                        }

                     /*
                     This block executes only for input non-vectorized
                     matrices with undefined column length. E.g. ckw01, ckw03
                     where inputs are (3,N) and (4,N). We uses the same check
                     as with vectorized inputs to ensure the required
                     consistency in column length.
                     */
                     if (extra.count == 0)
                        {
                        extra.first_vector_arg_index = ii;
                        extra.count                  = Dims[1];
                        }
                     else if (extra.count != Dims[1])
                        {
                        sprintf( msg,
                           "Input argument (`%s') must have "
                           "same length as `%s'",
                           argcheck[ii].name,
                           argcheck[extra.first_vector_arg_index].name);

                        mexErrMsgTxt(msg);
                        }

                     extra.vectorized[ii] = 1;

                     }
                  else
                     {

                     /*
                     Any other m, n value set signals an error.
                     */
                     sprintf( msg,
                        "MICE(BADARG): Input argument (`%s') with min_dims 2"
                        " improperly defined for this implementation of Mice."
                        " Error due to M=0 in an {M,N} ArgCheck assignment.",
                        argcheck[ii].name);

                     mexErrMsgTxt(msg);
                     }

                  break;

               default:

                  sprintf( msg,
                     "MICE(BUG): The min_dims value %d in argument `%s' not 0, "
                     "1, or 2. Check the argument specification in ArgCheck. "
                     "This indicates an interface bug. Contact NAIF.",
                      argcheck[ii].min_dims,
                      argcheck[ii].name);

                  mexErrMsgTxt(msg);
                  break;

               }

            }

         }
      else if (mxGetClassID(prhs[i]) == mxCHAR_CLASS)
         {

         /*
         Does this variable argument consist of string values?

         Vectorized string: N x M
         */
         if (argcheck[ii].min_dims != 0)
            {
            sprintf( msg,
                     "MICE(BUG): String input argument (`%s') "
                     "described with min_dims not zero. This "
                     "indicates an interface bug. Contact NAIF",
                     argcheck[ii].name);

            mexErrMsgTxt(msg);
            }

         /*
         In this case, we use the string length as the offset.
         Assume the variable as not vectorized till proven otherwise.
         */
         extra.offset[ii]     = Dims[1];
         extra.vectorized[ii] = 0;

         if( Dims[0] > 1 )
            {
            extra.vectorized[ii] = 1;
            }

         if (!argcheck[ii].is_vectorizable && extra.vectorized[ii])
            {
            sprintf( msg,
                     "MICE(BADARG): Input argument (`%s') "
                     "is not vectorizable.",
                     argcheck[ii].name);

            mexErrMsgTxt(msg);
            }

         if (argcheck[ii].is_vectorizable && extra.vectorized[ii])
            {

            /*
            Store the size of the input vector for return to the
            interface call.
            */
            if (extra.count == 0)
               {
               extra.first_vector_arg_index = ii;
               extra.count                  = Dims[0];
               }
            else if (extra.count != Dims[0] )
               {
               sprintf( msg,
                        "MICE(BADARG): Input argument (`%s') "
                        "must have the same vectorization measure as `%s'",
                        argcheck[ii].name,
                        argcheck[extra.first_vector_arg_index].name);

               mexErrMsgTxt(msg);
               }

            vec_flag = 1;

            }

         }
      else if (mxGetClassID(prhs[i]) == mxSTRUCT_CLASS)
         {

         /*
         Assume the variable as not vectorized till proven otherwise.

         Vectorized structure: 1 x N
         */
         extra.vectorized[ii] = 0;

         if( Dims[1] > 1 )
            {
            extra.vectorized[ii] = 1;
            }

         if (!argcheck[ii].is_vectorizable && extra.vectorized[ii])
            {
            sprintf( msg,
                     "MICE(BADARG): Input argument (`%s') "
                     "is not vectorizable.",
                     argcheck[ii].name);

            mexErrMsgTxt(msg);
            }

         if (argcheck[ii].is_vectorizable && extra.vectorized[ii])
            {

            /*
            Store the size of the input vector for return to the
            interface call.
            */
            if (extra.count == 0)
               {
               extra.first_vector_arg_index = ii;
               extra.count                  = Dims[1];
               }
            else if (extra.count != Dims[0] )
               {
               sprintf( msg,
                        "MICE(BADARG): Input argument (`%s') "
                        "must have the same vectorization measure as `%s'",
                        argcheck[ii].name,
                        argcheck[extra.first_vector_arg_index].name);

               mexErrMsgTxt(msg);
               }

            vec_flag = 1;

            }

         }
      else
         {

         sprintf( msg,
                  "MICE(BADARG): Input argument (`%s') "
                  "is not currently coded as an allowed input type.",
                  argcheck[ii].name);

         mexErrMsgTxt(msg);
         }

      /*
      Increment argcheck index.
      */
      ii++;
      }



   /*

   OUTPUT

   Check types, and allocate needed memory for output variables.
   Retrieve the dimension size array. Dims has length arg_dimension,
   each element of Dim containing the number of elements corresponding
   to that degree of freedom.
   */
   for ( i=0; i<nlhs; i++)
      {

      /*
      Output arg is of type double, int, or boolean. If the type flag has value
      MiceIgnore, do nothing with the output argument. The interface call will
      allocate the needed memory.

      n defines the number of parameters defined in 'sizearray',
      2 for scalars and vectors, 3 for vectorized matrices.

        - scalar [1,1]
        - vectorized scalar [1,R]
        - a vector [N,1]
        - vectorized vector [N,R]
        - matrix [N,N]
        - vectorized struct [1,R]
        - vectorized matrix [N,M,R], i.e. R versions of
          a [N,M] matrix.

      extra.count describes the measure of vectorization (R).
      */

      if ( argcheck[ii].min_dims != 0 &&
           argcheck[ii].min_dims != 1 &&
           argcheck[ii].min_dims != 2 )
         {
         sprintf( msg,
                  "MICE(BUG): The min_dims value in argument `%s' not 0, "
                  "1, or 2. Check the argument specification in ArgCheck. "
                  "This indicates an interface bug. Contact NAIF.",
                  argcheck[ii].name);

         mexErrMsgTxt(msg);
         }

      is_a_struct = argcheck[ii].type == MiceNameID   ||
                    argcheck[ii].type == MicePlane    ||
                    argcheck[ii].type == MiceEllipse  ||
                    argcheck[ii].type == MicePos      ||
                    argcheck[ii].type == MiceNear     ||
                    argcheck[ii].type == MiceSurf     ||
                    argcheck[ii].type == MicePool     ||
                    argcheck[ii].type == MiceSub_PS   ||
                    argcheck[ii].type == MiceSurf_PS  ||
                    argcheck[ii].type == MiceIlum     ||
                    argcheck[ii].type == MiceWnsumd   ||
                    argcheck[ii].type == MiceFrinfo   ||
                    argcheck[ii].type == MiceSFS      ||
                    argcheck[ii].type == MicePVN      ||
                    argcheck[ii].type == MiceState;

      if( argcheck[ii].type == MiceIgnore ||
          argcheck[ii].type == MiceChar    )
         {

         /*
         Do nothing with the output. The user will assign needed memory
         based on return parameters.
         */

         }
      else if (argcheck[ii].type == MiceDouble ||
               argcheck[ii].type == MiceInt    ||
               argcheck[ii].type == MiceBoolean)
         {

         if (argcheck[ii].is_vectorizable && vec_flag)
            {

            if (argcheck[ii].min_dims == 0)
               {
               n                = 2;
               sizearray[0]     = 1;
               sizearray[1]     = extra.count;

               extra.offset[ii] = 1;
               }
            else if (argcheck[ii].min_dims == 1)
               {
               n                = 2;
               sizearray[0]     = argcheck[ii].dims[0];
               sizearray[1]     = extra.count;

               extra.offset[ii] = argcheck[ii].dims[0];
               }
            else if (argcheck[ii].min_dims == 2)
               {
               n                = 3;
               sizearray[0]     = argcheck[ii].dims[0];
               sizearray[1]     = argcheck[ii].dims[1];
               sizearray[2]     = extra.count;

               extra.offset[ii] = argcheck[ii].dims[1] * argcheck[ii].dims[0];
               }

            extra.vectorized[ii] = 1;

            }
         else
            {

            switch ( argcheck[ii].min_dims )
               {
               case 0:
                  sizearray[0] = 1;
                  sizearray[1] = 1;
                  break;

               case 1:
                  sizearray[0] = argcheck[ii].dims[0];
                  sizearray[1] = 1;
                  break;

               case 2:
                  sizearray[0] = argcheck[ii].dims[0];
                  sizearray[1] = argcheck[ii].dims[1];
                  break;

               default:
                  sprintf(msg,
                     "MICE(BUG): The min_dims value in argument `%s' not 0, "
                     "1, or 2. Check the argument specification in ArgCheck. "
                     "This indicates an interface bug. Contact NAIF.",
                      argcheck[ii].name);

                  mexErrMsgTxt(msg);
                  break;
               }

            n = 2;
            extra.vectorized[ii] = 0;

            }

         /*
         Allocate output memory of the proper size and type for dps, ints,
         and bools.
         */
         if (argcheck[ii].type == MiceDouble)
            {
            plhs[i] = mxCreateNumericArray( n,
                                            sizearray,
                                            mxDOUBLE_CLASS,
                                            mxREAL);
            }
         else if ( argcheck[ii].type == MiceInt    ||
                   argcheck[ii].type == MiceBoolean)
            {
            plhs[i] = mxCreateNumericArray( n,
                                            sizearray,
                                            mxINT32_CLASS,
                                            mxREAL);
            }
         else
            {

            /*
            This should never signal.
            */
            mexErrMsgTxt( "MICE(BUG): Serious error! Output type ID error. "
                          "This indicates an interface bug. Contact NAIF.");
            }

         }

      else if ( is_a_struct )
         {

         if (argcheck[ii].is_vectorizable && vec_flag)
            {

            n                = extra.count;
            extra.offset[ii] = 1;

            }
         else
            {
            n = 1;
            }

         nfields = 0;
         fnames  = NULL;
         ftypes  = NULL;
         fsizes  = NULL;

         struct_fields( argcheck[ii].type,
                                    &nfields, &fnames, &ftypes, &fsizes);

         /*
         Create mxArray for the struct.
         */
         plhs[i] = mxCreateStructMatrix(1, n, nfields, fnames);

         for (j=0; j<n; j++)
            {

            /*
            Create mxArrays for each of the fields. Note, memory
            allocation for each member of the struct (important to
            know this).
            */
            for (k=0;k<nfields;k++)
               {

               if ( ftypes[k]==MiceDouble || ftypes[k]==MiceInt )
                  {
                  mxSetField( plhs[i],
                              j,
                              fnames[k],
                              mxCreateDoubleMatrix(fsizes[k], 1, mxREAL));
                  }
               else if (ftypes[k]==MiceBoolean )
                  {
                  mxSetField(plhs[i],
                             j,
                             fnames[k],
                             mxCreateLogicalMatrix(fsizes[k], 1) );
                  }
               else if ( ftypes[k]==MiceChar )
                  {
                  sizearray[0] = 1;
                  sizearray[1] = DEFAULT_STR_LENGTH;

                  mxSetField(plhs[i],
                             j,
                             fnames[k],
                             mxCreateCharArray(2, sizearray));
                  }

               }

            }

         }
      else
         {
         sprintf( msg,
                  "MICE(BUG): Return argument (`%s') of unknown type. "
                  "This indicates an interface bug. Contact NAIF.",
                  argcheck[ii].name);

         mexErrMsgTxt(msg);
         }

      /* Increment argcheck index. */
      ii++;
      }

   return &extra;
   }





