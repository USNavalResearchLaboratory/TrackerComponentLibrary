%-Abstract
%
%   MiceOccult.m declares occultation specific definitions. 
%   
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
% 
%
%-I/O
%
%   None.
% 
%   the call:
% 
%       MiceOccult
%
%-Examples
%
%   None.
% 
%-Particulars
%
%    The following integer codes indicate the geometric relationship
%    of the three bodies.
% 
%    The meaning of the sign of each code is given below.
% 
%        Code sign     Meaning
%        ---------     ------------------------------
%           > 0        The second ellipsoid is
%                      partially or fully occulted
%                      by the first.
% 
%           < 0        The first ellipsoid is
%                      partially of fully
%                      occulted by the second.
% 
%           = 0        No occultation.
% 
%    The meanings of the codes are given below. The variable names indicate
%    the type of occultation and which target is in the back. For example, 
%    SPICE_OCCULT_TOTAL1 represents a total occultation in which the 
%    first target is in the back (or occulted by) the second target.
% 
%        Name                       Code     Meaning
%        ------                     -----    ------------------------------
%        SPICE_OCCULT_TOTAL1         -3      Total occultation of first
%                                            target by second. First target
%                                            is in back.
% 
%        SPICE_OCCULT_ANNLR1         -2      Annular occultation of first
%                                            target by second.  The second
%                                            target does not block the limb is
%                                            in of the first. First target
%                                            back.
% 
%        SPICE_OCCULT_PARTL1         -1      Partial occultation of first
%                                            target by second target. First 
%                                            target is in back.
% 
%        SPICE_OCCULT_NOOCC           0      No occultation or transit: both
%                                            objects are completely visible
%                                            to the observer.
% 
%        SPICE_OCCULT_PARTL2          1      Partial occultation of second
%                                            target by first target. Second 
%                                            target is in back.
% 
%        SPICE_OCCULT_ANNLR2          2      Annular occultation of second
%                                            target by first. Second target is
%                                            in back.
% 
%        SPICE_OCCULT_TOTAL2          3      Total occultation of second
%                                            target by first. Second target is
%                                            in back.
%
%-Required Reading
%
%   None.
%
%-Version
%
%   -Mice Version 1.0.0, 20-FEB-2014, SCK (JPL)
%
%-Index_Entries
%
%   Include occultation code parameters
%
%-&


global SPICE_OCCULT_TOTAL1  ...
       SPICE_OCCULT_ANNLR1  ...
       SPICE_OCCULT_PARTL1  ...
       SPICE_OCCULT_NOOCC   ...
       SPICE_OCCULT_PARTL2  ...
       SPICE_OCCULT_ANNLR2  ...
       SPICE_OCCULT_TOTAL2


SPICE_OCCULT_TOTAL1 = -3;
SPICE_OCCULT_ANNLR1 = -2;
SPICE_OCCULT_PARTL1 = -1;
SPICE_OCCULT_NOOCC  =  0;
SPICE_OCCULT_PARTL2 =  1;
SPICE_OCCULT_ANNLR2 =  2;
SPICE_OCCULT_TOTAL2 =  3;

