The library has been modified by making line_search_backtracking,
line_search_backtracking_owlqn, and line_search_morethuente non-static and
moved prototypes into lbfgs.h so the line search methods can be used by
other functions. Also, the declaration of tag_callback_data into lbfgs.h
and checks for MATLAB_MEX_FILE were added to the code dealing with memory
allocation and deallocation to ease the integration into Matlab.

Additionally, some files that are not used have been removed. These are
primarily files that were subject to a different license with undesirable
copyleft clauses. The retained library is copyleft-free.

The library is originally from
http://www.chokkan.org/software/liblbfgs/

July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.