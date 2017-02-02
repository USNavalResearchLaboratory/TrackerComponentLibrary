/*
-Procedure

   cspice_params.h

-Abstract

    Parameter values from SPICELIB code for use in CSPICE.

    KEEP THE ASSIGNMENTS DEFINED IN THIS FILE SYNCHED WITH 
    THE CORRESPONDING ASSIGNMENTS IN POOL.F, ZZRVAR.F, DAFFA.F, 
    ERRHND.INC.

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

   Paramters

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

   -Paramter values taken from pool.f version:

       SPICELIB Version 10.1.0, 20-JUN-2013 (BVS)

         MAXVAR      is the maximum number of variables that the
                     kernel pool may contain at any one time.
                     MAXVAR should be a prime number.

         MAXLEN      is the maximum length of the variable names
                     that can be stored in the kernel pool.

         MAXVAL      is the maximum number of distinct values that
                     may belong to the variables in the kernel pool.
                     Each variable must have at least one value, and
                     may have any number, so long as the total number
                     does not exceed MAXVAL. MAXVAL must be at least
                     as large as MAXVAR.

         MXNOTE      is the maximum number of distinct variable-agents
                     pairs that can be maintained by the kernel pool.
                     (A variable is "paired" with an agent, if that agent
                     is to be notified whenever the variable is updated.)

         MAXAGT      is the maximum number of agents that can be kept
                     on the distribution list for notification of updates
                     to kernel variables.

         MAXCHR      is the maximum number of characters that can be
                     stored in a component of a string valued kernel
                     variable.

         MAXLIN      is the maximum number of character strings that
                     can be stored as data for kernel pool variables.

   -Paramter values taken from zzrvar.f version:

      SPICELIB Version 1.6.0, 06-AUG-2002 (BVS)

         LINLEN      is the maximum length of a line in the kernel file.

   -Paramter values taken from daffa.f version:
    
      SPICELIB Version 3.0.0, 16-NOV-2001 (FST)

         MAXNDC

         MAXNIC

   -Paramter values taken from dafah.f version:

      SPICELIB Version 9.0.0, 09-NOV-2006 (NJB)

         MAXSUM

      Paramter SIDLEN derived from:

                             (NI + 1)
         SIDLEN = 8 * ( ND + -------- )     (Note that this is
                                2           integer division.)

          with ND = 2, and NI = 6 for SPKs.

   -Parameter values take from errhnd.h version:
   
      SPICELIB Version 3.0.0, 14-JAN-2013 (EDW)

         LMSGLN

         SMSGLN

-Examples

    None.

-Restrictions

    None.

-Literature_References

    DAF.REQ
    ERROR.REQ
    POOL.REQ

-Author_and_Institution

    E. D. Wright    (JPL)

-Version

   -CSPICE Version 1.3.0, 18-NOV-2013 (EDW) (BVS)

      Added parameters:
      
         MAXNDC
         MAXNIC
         MAXSUM
         SIDLEN
         LMSGLN
         SMSGLN

      Updated to parameter values to match those defined
      in pool.f:
      
         MAXVAR to 26003
         MAXVAL to 400000
         MAXLIN to 15000
         MXNOTE to (MAXVAR * 5)

   -CSPICE Version 1.2.0, 24-MAY-2010 (EDW) (NJB)

      Increased MAXVAL to 200000.
        
   -CSPICE Version 1.1.0, 23-FEB-2009 (EDW)

      Added LINLEN parameter.

   -CSPICE Version 1.0.0, 27-APR-2006 (EDW)

-Index_Entries

   parameter definitions

-&
*/

#define         MAXVAR         26003

#define         MAXVAL         400000

#define         MAXLIN         15000

#define         MAXCHR         80

#define         MXNOTE         (MAXVAR * 5)

#define         MAXLEN         32
 
#define         MAXAGT         1000 

#define         LINLEN         132

#define         MAXNDC         124

#define         MAXNIC         250

#define         MAXSUM         125

#define         SIDLEN         40

#define         LMSGLN         23 * 80
   
#define         SMSGLN         25
