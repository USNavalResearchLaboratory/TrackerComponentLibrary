/**UTC2DAYCOUNTC A C-function for converting a universal coordinated time
 *               (UTC) given as a two-part pseudo-Julian to a date in terms
 *               of the year, day number and second in the day.
 *
 *INPUTS:    Jul1, Jul2  Two parts of a Julian date given in UTC. The units
 *                       of the date are days. The full date is the sum of
 *                       both terms. The date is broken into two parts to
 *                       provide more bits of precision. It does not matter
 *                       how the date is split.
 *
 *          year    A pointer to an integer to hold the year in the
 *                  Gregorian calendar under UTC time.
 *         dayCount A pointer to an integer to hold the day of the year,
 *                  accounting for leap years. The returned value will be
 *                  >=1.
 *          second  A pointer to a double to hold the second of the day,
 *                  accounting for leap seconds. The returned value will be
 *                  >=0.
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
 *Leap years are taken into account. February has 29 days if the year is
 *divisible by 4, unless the year is divisible by 100 unless the year is
 *divisible by 400.
 *
 *UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
 *be called with an earlier date.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "SharedTimeCode.h"

int UTC2DayCountC(double Jul1, double Jul2, int *year, int *dayCount, double *second) {
    int retVal,month,day;
    double dayFrac,ddt;

    retVal=iauJd2cal(Jul1, Jul2, year, &month, &day, &dayFrac);
    
    if(retVal!=0&&retVal!=1) {
        return retVal;//A date conversion error has occurred.
    }

    //Compute the day number in the year accounting for whether it is a
    //leap year according to the Gregorian calendar.
    *dayCount=day;
    switch(month) {
        case 12:
            *dayCount+=30;
        case 11:
            *dayCount+=31;
        case 10:
            *dayCount+=30;
        case 9:
            *dayCount+=31;
        case 8:
            *dayCount+=31;
        case 7:
            *dayCount+=30;
        case 6:
            *dayCount+=31;
        case 5:
            *dayCount+=30;
        case 4:
            *dayCount+=31;
        case 3:
            //The rule for leap year is: a year that is divisible by 4
            //is a leap year unless, the year is divisible by 100
            //unless that year is divisible by 400.
            if(*year%4==0) {
                if(*year%100==0&&*year%400!=0) {
                    *dayCount+=28;//Not a leap year.
                }
                else {
                    *dayCount+=29;//It is a leap year.
                }
            } else {//Not a leap year.
                *dayCount+=28;
            }
        case 2:
            *dayCount+=31;
        default:
            break;
    }
    
//Now, figure out the second in the day accounting for leap seconds.
    
/*The following few lines are based on code in dtf2d.c in the IAU's library
 *to determine whether the day in question contains a leap second.*/
    {
      double dat1, dat2, w;
      int nextYear,nextMonth,nextDay;
   /* TAI-UTC today. */
      iauDat(*year, month, day, 0.0, &dat1);
   /* TAI-UTC tomorrow. */
      iauJd2cal (Jul1+Jul2, 1.0, &nextYear, &nextMonth, &nextDay, &w);
      iauDat(nextYear, nextMonth, nextDay, 0.0, &dat2);
   /* The change in TAI-UTC (seconds). */
      ddt = dat2 - dat1;
    }
      
/*The final minute in the day contains 60+ddt seconds.*/
    *second=dayFrac*(23*60*60+59*60+(60+ddt));
    
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
