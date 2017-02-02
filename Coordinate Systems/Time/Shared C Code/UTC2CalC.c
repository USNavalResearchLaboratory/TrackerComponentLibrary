/**UTC2CALC A C-function for converting a universal coordinated time (UTC)
 *          given as a two-part pseudo-Julian to a date in terms of the
 *          Gregorian calendar in years, months, days, hours, minutes and
 *          seconds. 
 *
 *INPUTS:    Jul1, Jul2  Two parts of a Julian date given in UTC. The units
 *                       of the date are days. The full date is the sum of
 *                       both terms. The date is broken into two parts to
 *                       provide more bits of precision. It does not matter
 *                       how the date is split.
 *
 *          year    A pointer to an integer to hold the year in the
 *                  Gregorian calendar under UTC time.
 *          month   A pointer to an integer to hold the month in the
 *                  Gregorian calendar under UTC time.
 *          day     A pointer to an integer to hold the day in the
 *                  Gregorian calendar under UTC time.
 *          hour    A pointer to an integer to hold the hour the Gregorian
 *                  calendar under UTC time.
 *          minute  A pointer to an integer to hold the minute in the
 *                  Gregorian calendar under UTC time.
 *          second  A pointer to a double to hold the second in the
 *                  Gregorian calendar under UTC time. This includes the
 *                  possibility of a leap second on the day in question.
 *
 *OUTPUTS:  The outputs are placed in the passed pointers for the date. The
 *          return value is that of the iauJd2cal function in the SOFA
 *          library. Thus, 0 is fine, 1 is a dubious date was entered and
 *          -1 means that an unacceptable date was entered.
 *
 *This makes use of the function iauJd2cal in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library
 *whereby the output from that function is transformed to be in terms of
 *hours minutes and seconds rather than as a fraction of a day.
 *The day from which the fraction is taken is assumed to have 24 hours,
 *with 60 minutes per hour and 60 seconds per minute.
 *
 *The UTC date is only pseudo-Julian, because there is not a fixed number
 *of seconds in a Julian day. The convention used in the IAU standard is
 *that the Julian day matches the UTC day regardless of whether the UTC day
 *is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
 *seconds).
 *
 *UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
 *be called with an earlier date.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "SharedTimeCode.h"

int UTC2CalC(double Jul1, double Jul2,int *year, int *month, int *day, int *hourRet, int *minuteRet, double *second) {
    double dayFrac,ddt, hour, minute;
    int retVal;
    
    /*Call the function in the SOFA library.*/
    retVal=iauJd2cal(Jul1, Jul2, year, month, day, &dayFrac);
    
    if(retVal!=0&&retVal!=1) {
        return retVal;
    }

/*The following few lines are based on code in dtf2d.c in the IAU's library
 *to determine whether the day in question contains a leap second.*/
    {
      double dat1, dat2, w;
      int nextYear,nextMonth,nextDay;
   /* TAI-UTC today. */
      iauDat(*year, *month, *day, 0.0, &dat1);
   /* TAI-UTC tomorrow. */
      iauJd2cal (Jul1+Jul2, 1.0, &nextYear, &nextMonth, &nextDay, &w);
      iauDat(nextYear, nextMonth, nextDay, 0.0, &dat2);
   /* The change in TAI-UTC (seconds). */
      ddt = dat2 - dat1;
    }
      
/*The final minute in the day contains 60+ddt seconds. Given the number of
 *seconds in a day, back out the number of hours and minutes.*/
    *second=dayFrac*(23*60*60+59*60+(60+ddt));
    {
        double temp=floor(*second/(60*60));
        
        if(temp>23) {
            hour=23;
        } else {
            hour=temp;
        }
    }
    
    *second=*second-hour*60*60;
    {
        double temp=floor(*second/60);
        if(temp>59) {
            minute=59;
        } else {
            minute=temp;
        }
    }
    *second=*second-minute*60;
    *hourRet=(int)hour;
    *minuteRet=(int)minute;
    
    return retVal;
}

/*LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
