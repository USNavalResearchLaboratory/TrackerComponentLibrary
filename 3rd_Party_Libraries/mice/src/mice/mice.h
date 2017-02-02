/*

-Abstract

   Define and assign the Mice particular variable types.
   
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

#ifndef mice_h
#define mice_h

/*
   Define needed globals.

   1) Size for temporary strings in interface calls.

   2) SCALAR has value -2 to compensate for the difference between array index
      bases in MATLAB (base 1) and C (base 0). The zzerror routine determines 
      whether to return a vector based error string or a scalar base string 
      based on the value of 'cnt' as less than zero. The expression 
      'cnt + INDEX_BASE' passes to zzerror, INDEX_BASE having value either
      0 or 1. If 1, SCALAR requires a value less than -1 to satisfy zzerror 
      logic.

*/
   
#define DEFAULT_STR_LENGTH   1024
#define SCALAR               -2

/*
   Note, MATLAB uses base 1 indexing (as is proper), while C uses base 0
   (corresponding to an offset). Add 1 to 'cnt'.
*/
   
#define INDEX_BASE          1L

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


   -Mice Version 1.0.2 18-JUN-2012 (EDW)

      Removed unused macro S_DBL_RET_ARGV.

   -Mice Version 1.0.1 23-FEB-2009 (EDW)

      Added the *_OFF, *_IN, *_OUT, index aliases.

   -Mice Version 1.0.0 31-JAN-2008 (EDW)

-Index_Entries

   MATLAB

-&
*/



/*
Prototypes.
*/
mxArray * zzmice_CreateIntScalar( SpiceInt n);


/*
Variable's macros.
*/
#define       S_INT_ARGV(x)         *(SpiceInt*)mxGetData(prhs[x])

#define       A_INT_ARGV(x)          (SpiceInt*)mxGetData(prhs[x])

#define       S_DBL_ARGV(x)         *(SpiceDouble*)mxGetData(prhs[x])

#define       A_DBL_ARGV(x)          (SpiceDouble*)mxGetData(prhs[x])

#define       A_DBL_RET_ARGV(x)      (SpiceDouble*)mxGetData(plhs[x])

#define       A_BOOL_RET_ARGV(x)     (SpiceBoolean*)mxGetData(plhs[x])

#define       A_INT_RET_ARGV(x)      (SpiceInt*)mxGetData(plhs[x])


/*
Index aliases.
*/
#define       ONE_OFF            0
#define       TWO_OFF            1
#define       THREE_OFF          2
#define       FOUR_OFF           3
#define       FIVE_OFF           4
#define       SIX_OFF            5
#define       SEVEN_OFF          6
#define       EIGHT_OFF          7

#define       ONE_IN             1
#define       TWO_IN             2
#define       THREE_IN           3
#define       FOUR_IN            4
#define       FIVE_IN            5
#define       SIX_IN             6
#define       SEVEN_IN           7
#define       EIGHT_IN           8

#define       ONE_OUT            0
#define       TWO_OUT            1
#define       THREE_OUT          2
#define       FOUR_OUT           3
#define       FIVE_OUT           4
#define       SIX_OUT            5
#define       SEVEN_OUT          6
#define       EIGHT_OUT          7



/*
Simple macro based on ALLOC_CHECK to ensure a zero value alloc count
at end of routine, if not, pass the error message to the MATLAB interpreter.
Note, the need to use this macro exists only in those routines
allocating/deallocating memory.
*/
#define MICE_ALLOC_CHECK                                                   \
             if ( alloc_count() != 0 )                                     \
                {                                                          \
                setmsg_c ( "MICE(BUG): Malloc/Free count not zero at end " \
                           "of routine. Malloc count = #. Contact NAIF." );\
                errint_c ( "#", alloc_count ()    );                       \
                sigerr_c ( "SPICE(MALLOCCOUNT)"        );                  \
                mice_fail(SCALAR);                                         \
                }



/*
This checks failed_c to detect a SPICE error signal.

Code should call this macro after a call to any subroutine or function
that might signal a SPICE error.
*/
#define CHECK_CALL_FAILURE(cnt)  \
             if ( failed_c())    \
                {                \
                mice_fail(cnt);  \
                }


/*
CHECK_CALL_FAILURE_MEM(n,arr) performs the same functions as in
CHECK_CALL_FAILURE but also frees memory blocks allocated by the
alloc_SpiceString_C_array routines.
*/
#define CHECK_CALL_FAILURE_MEM(n,arr)                \
             if ( failed_c())                        \
                {                                    \
                free_SpiceString_C_array ( n, arr ); \
                mice_fail(SCALAR);                   \
                }


/*
CHECK_CALL_FAILURE_MEM1(cnt,n,arr,pntr) performs the same functions as in
CHECK_CALL_FAILURE but also frees memory blocks allocated by the
alloc_SpiceString_C_array/alloc_SpiceString_Pointer_array routines.
*/
#define CHECK_CALL_FAILURE_MEM1(cnt,n,arr,pntr)      \
             if ( failed_c())                        \
                {                                    \
                free_SpiceString_C_array ( n, arr ); \
                free_SpiceMemory         ( pntr   ); \
                mice_fail(cnt);                      \
                }


#endif

