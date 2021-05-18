/* This is a suite of C and C++ functions for type checking inputs to mex
 * files and performing other helper actions when interfacing Matlab and C
 * or C++ code.
 *
 *The following functions that can be used in both C and C++ programs:
 *checkRealArray            Gives an error if an array is not real, has
 *                          more than two indices, or is empty.
 *checkRealDoubleArray      Gives an error if an array is not composed of
 *                          real doubles, has more than two indices, or is
 *                          empty.
 *checkDoubleArray          Gives an error if an array is not composed of
 *                          doubles, has more than two indices, or is
 *                          empty.
 *checkRealDoubleHypermatrix This is the same as checkRealDoubleArray
 *                          except it does not throw an error if the matrix
 *                          has more than two indices.
 *checkRealInt32Array       Gives an error if an array is not composed of
 *                          real int32_t values, has more than two indices,
 *                          or is empty.
 *checkInt32Array           Gives an error if an array is not composed of
 *                          int32_t values, has more than two indices, or
 *                          is empty.
 *verifySizeReal            Verify that a Matlab matrix is 2D, has the
 *                          desired dimensions, and is real.
 *convert2DReal2DoubleMat   Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          doubles. The returned Matrix should be freed
 *                          with mxDestroyArray. the input matrix is not
 *                          modified.
 *convert2DReal2SignedIntMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          ints. The returned Matrix should be freed
 *                          with mxDestroyArray. the input matrix is not
 *                          modified.
 *convert2DReal2UnsignedIntMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          unsigned ints. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *convert2DReal2SignedSizeMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          ptrdiff_t. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *convert2DReal2UnsignedSizeMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          size_t. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *getIntFromMatlab          Turns a passed Matlab data type into a real
 *                          integer. Provides an error if the Matlab data
 *                          is of a type that can not be transformed into
 *                          an integer or if it is complex.
 *getSizeTFromMatlab        Turns a passed Matlab data type into a
 *                          positive, real size_t integer value. Provides
 *                          an error if the Matlab data is of a type that
 *                          can not be transformed into an size_t value or
 *                          if it is negative or complex.
 *getPtrDiffTFromMatlab     Turns a passed Matlab data type into a real
 *                          ptrdiff_t integer value. Provides
 *                          an error if the Matlab data is of a type that
 *                          can not be transformed into an ptrdiff_t value.
 *copySizeTArrayFromMatlab  Returns a copy of an array from Matlab as an
 *                          array of size_t values and indicates the length
 *                          of the array by modifying an input parameter.
 *                          The array is allocated using mxAlloc and thus
 *                          should be freed using mxFree.
 *copyPtrDiffTArrayFromMatlab  Returns a copy of an array from Matlab as an
 *                          array of ptrdiff_t values and indicates the
 *                          length of the array by modifying an input
 *                          parameter.  The array is allocated using
 *                          mxAlloc and thus should be freed using mxFree.
 *copyBoolArrayFromMatlab   Returns a copy of an array from Matlab
 *                          converted into boolean variables and indicates
 *                          the length of the array by modifying an input
 *                          parameter. The array is allocated using mxAlloc
 *                          and thus should be freed using mxFree.
 *getDoubleFromMatlab       Turns a passed Matlab data type into a real
 *                          double variable. Provides an error if the
 *                          Matlab data is of a type that can not be
 *                          transformed into a double or if it is
 *                          complex.
 *getBoolFromMatlab         Turns a passed Matlab data type into a real
 *                          boolean variable. Provides an error if the
 *                          Matlab data is of a type that can not be
 *                          transformed into a boolean or if it is
 *                          complex.
 *getInt32FromMatlab        Turns a passed Matlab data type into an real
 *                          int32_t variable. Provides an error if the
 *                          Matlab data is of a type that can not be
 *                          transformed into an int32_t or if it is
 *                          complex.
 *doubleMat2Matlab          Convert an array or matrix of doubles into
 *                          a matrix to be returned to Matlab.
 *floatMat2MatlabDoubles    Convert an array or matrix of single-precision
 *                          floats into  a matrix to be returned to Matlab
 *                          as double floating point numbers (the default
 *                          data format in Matlab).
 *intMat2MatlabDoubles      Convert an array or matrix of integers into
 *                          a matrix to be returned to Matlab as double
 *                          floating point numbers (the default data format
 *                          in Matlab).
 *int32Mat2MatlabDoubles    Convert an array or matrix of int32_T values
 *                          into a matrix to be returned to Matlab as
 *                          double floating point numbers (the default data
 *                          format in Matlab).
 *sizeTMat2MatlabDoubles    Convert an array or matrix of size_t values
 *                          into a matrix to be returned to Matlab as
 *                          double floating point numbers (the default data
 *                          format in Matlab).
 *ptrDiffTMat2MatlabDoubles Convert an array or matrix of ptrdiff_t values
 *                          into a matrix to be returned to Matlab as
 *                          double floating point numbers (the default data
 *                          format in Matlab).
 *boolMat2Matlab            Convert an array or matrix of bool values into
 *                          a matrix to be returned to Matlab (as type
 *                          mxLogical).
 *allocUnsignedSizeMatInMatlab Allocate a matrix for holding size_t values
 *                          that can be used in Matlab.
 *allocSignedSizeMatInMatlab Allocate a matrix for holding ptrdiff_t values
 *                          that can be used in Matlab.
 *allocUnsignedCharMatInMatlab Allocate a matrix for holding unsigned char
 *                          values that can be used in Matlab.
 *allocUnsignedIntMatInMatlab Allocate a matrix for holding unsigned int
 *                          values that can be used in Matlab.
 *allocSignedIntMatInMatlab Allocate a matrix for holding int values
 *                          that can be used in Matlab.
 *unsignedSizeMat2Matlab    Convert an array or matrix of size_t values
 *                          into a matrix to be returned by Matlab.
 *unsignedCharMat2Matlab    Convert an array or matrix of unsigned char
 *                          values into a matrix to be returned by Matlab.
 *signedSizeMat2Matlab      Convert an array or matrix of ptrdiff_t values
 *                          into a matrix to be returned to Matlab.
 *getScalarMatlabClassConst Get a double scalar constant that is defined in
 *                          a Matlab class.
 *pointerIsAligned          Determine whether a pointer is aligned to a
 *                          certain number of bytes. This is useful if one
 *                          wishes to use SIMD instructions in the
 *                          processor,because Matlab does not specify that
 *                          matrices will maintain a 16-byte (or other)
 *                          alignment, which is generally necessary for
 *                          such instructions to be particularly fast.
 *
 *The following functions are for use when interfacing C and C++ programs:
 *ptr2Matlab                Convert a pointer to a class instance into a
 *                          data type that can be returned to Matlab so
 *                          that the class can be used during another call
 *                          to a mex function.
 *Matlab2Ptr                Convert a Matlab value obtained from ptr2Matlab
 *                          back to a pointer to a C++ class.
 *
 *The following functions are only available in C++ as they use template
 *classes:
 *mat2MatlabDoubles         Convert an array or matrix of some type to a
 *                          matrix of doubles to be returned in Matlab. The
 *                          components of the array must be convertible to
 *                          doubles using the static_cast function. The
 *                          use is the same as the type-specific functions
 *                          such as intMat2MatlabDoubles.
 *getComplexDoubleFromMatlab Given a scalar numeric data types in Matlab,
 *                          return a complex<double> with the real and
 *                          imaginary parts of the input. If the input was
 *                          purely imaginary or real, then the
 *                          corresponding part of the return value will be
 *                          zero. This function provides an error if the
 *                          data type is not a type that can be converted
 *                          to a double or is not scalar.
 *
 *September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MEXHELP
#define MEXHELP
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
//For memcpy
#include <cstring>
//Defines the size_t and ptrdiff_t types
#include <cstddef>
#include <cstdint>
//For functions to import comples values.
#include <complex>
#else
//For memcpy
#include <string.h>
//Defines the size_t and ptrdiff_t types
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

//This is needed for the bool type to be defined in C.
//C99 has stdbool, earlier versions do not.
#if __STDC_VERSION__>=199901L

//If using *NIX or Mac OS X, the functions should be defined without
//anything extra in math.h
#include<math.h>
#define fMax(a,b) fmax(a,b)
#define fMin(a,b) fmin(a,b)
#define isFinite(x) isfinite(x)
#define copySign(a,b) copysign(a,b)

#else

//These functions are defined in the C standard, but do not seem to be
//supported by Microsoft's compiler, so they are defined here with
//different capitalization.
double fMax(double a,double b) {
//Defined in the C99 standard, but not supported by Microsoft.
    if(a>b) {
        return a;
    }
    else {
        return b;
    }
}

double fMin(double a,double b) {
//Defined in the C99 standard, but not supported by Microsoft.
    if(a>b) {
        return b;
    }
    else {
        return a;
    }
}

//Microsoft does not declare isfinite and copysign (in the C99 standard),
//but it does declare _finite and _copysign in float.h for some
//inexplicable reason.
#include <float.h>
#define isFinite(x)_finite(x)
#define copySign(a,b)_copysign(a,b)

#endif
#endif
//void checkRealArray(const mxArray * const val);
//void checkRealDoubleArray(const mxArray * const val);
//void checkDoubleArray(const mxArray * const val);
//void checkRealDoubleHypermatrix(const mxArray * const val);
//void checkRealInt32Array(const mxArray * const val);
//void checkInt32Array(const mxArray * const val);
//void verifySizeReal(const size_t M, const size_t N, const mxArray * const val);
//mxArray *convert2DReal2DoubleMat(const mxArray * const val);
//mxArray *convert2DReal2SignedIntMat(const mxArray * const val);
//mxArray *convert2DReal2UnsignedIntMat(const mxArray * const val);
//mxArray *convert2DReal2SignedSizeMat(const mxArray * const val);
//mxArray *convert2DReal2UnsignedSizeMat(const mxArray * const val);

//int getIntFromMatlab(const mxArray * const val);
//size_t getSizeTFromMatlab(const mxArray * const val);
//ptrdiff_t getPtrDiffTFromMatlab(const mxArray * const val);
//size_t *copySizeTArrayFromMatlab(const mxArray * const val, size_t *arrayLen);
//ptrdiff_t *copyPtrDiffTArrayFromMatlab(const mxArray * const val, size_t *arrayLen);
//bool *copyBoolArrayFromMatlab(const mxArray * const val, size_t *arrayLen);
//double getDoubleFromMatlab(const mxArray * const val);
//bool getBoolFromMatlab(const mxArray * const val);
//int32_t getInt32FromMatlab(const mxArray * const val);

//mxArray *doubleMat2Matlab(const double * const arr,const size_t numRow, const size_t numCol);
//mxArray *floatMat2MatlabDoubles(const float * const arr,const size_t numRow, const size_t numCol);
//mxArray *intMat2MatlabDoubles(const int * const arr,const size_t numRow, const size_t numCol);
//mxArray *sizeTMat2MatlabDoubles(const size_t * const arr,const size_t numRow, const size_t numCol);
//mxArray *ptrDiffTMat2MatlabDoubles(const ptrdiff_t * const arr,const size_t numRow, const size_t numCol);
//mxArray *boolMat2Matlab(const bool * const arr,const size_t numRow, const size_t numCol);
//mxArray *allocUnsignedSizeMatInMatlab(const size_t numRow, const size_t numCol);
//mxArray *allocSignedSizeMatInMatlab(const size_t numRow, const size_t numCol);
//mxArray *allocUnsignedCharMatInMatlab(const size_t numRow, const size_t numCol);
//mxArray *allocUnsignedIntMatInMatlab(const size_t numRow, const size_t numCol);
//mxArray *allocSignedIntMatInMatlab(const size_t numRow, const size_t numCol);
//mxArray *unsignedSizeMat2Matlab(const size_t * const arr, const size_t numRow, const size_t numCol);
//mxArray *unsignedCharMat2Matlab(const unsigned char * const arr, const size_t numRow, const size_t numCol);
//mxArray *signedSizeMat2Matlab(const ptrdiff_t * const arr, const size_t numRow, const size_t numCol);
//double getScalarMatlabClassConst(const char * const className, const char * const constName);
//bool pointerIsAligned(const void *pointer, const size_t byteCount);

#ifdef __cplusplus
//One cannot use prototypes for template functions. The template C++-only
//functions are full implementations.
std::complex <double> getComplexDoubleFromMatlab(const mxArray * const val);
#endif

void checkRealArray(const mxArray * const val) {
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
}

void checkDoubleArray(const mxArray * const val) {
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    if(mxGetClassID(val)!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("A parameter that should be a real double is of a different data type.");
    }
}

void checkRealDoubleArray(const mxArray * const val){
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    checkDoubleArray(val);
}

void checkRealDoubleHypermatrix(const mxArray * const val){
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)<1){
        mexErrMsgTxt("A parameter has too few dimensions.");
    }
    
    if(mxGetClassID(val)!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("A parameter that should be a real double is of a different data type.");
    }
}

void checkInt32Array(const mxArray * const val) {
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of int32_t values is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    if(mxGetClassID(val)!=mxINT32_CLASS) {
        mexErrMsgTxt("A parameter that should be int32_t values is of a different data type.");
    }
}

void checkRealInt32Array(const mxArray * const val) {
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of int32_t values has complex components.");
    }
    
    checkInt32Array(val);
}

void verifySizeReal(const size_t M, const size_t N, const mxArray * const val) {
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix has complex components.");
    }
    
    if(mxGetM(val)!=M||mxGetN(val)!=N) {
        mexErrMsgTxt("A parameter has the wrong dimensionality.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
}

mxArray *convert2DReal2DoubleMat(const mxArray * const val) {
    mxArray *retMat;
    double *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    retMat=mxCreateDoubleMatrix(M,N,mxREAL);
    retData=mxGetDoubles(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *data=mxGetChars(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *data=mxGetLogicals(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=mxGetDoubles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=mxGetSingles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=mxGetInt8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=mxGetUint8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=mxGetInt16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=mxGetUint16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=mxGetInt32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=mxGetUint32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=mxGetInt64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=mxGetUint64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that cannot be used.");
    }
    
    return retMat;
}

mxArray *allocSignedIntMatInMatlab(const size_t numRow, const size_t numCol){
    mxArray *retMat=NULL;
    
    switch(sizeof(int)) {
        case 4:
            retMat=mxCreateNumericMatrix(numRow,numCol,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat=mxCreateNumericMatrix(numRow,numCol,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *convert2DReal2SignedIntMat(const mxArray * const val) {
    mxArray *retMat;
    int *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocSignedIntMatInMatlab(M,N);
    if(sizeof(int)==4) {//32 bit
        retData=(int*)mxGetInt32s(retMat);
    } else {//64 bit
        retData=(int*)mxGetInt64s(retMat);
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *data=mxGetChars(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *data=mxGetLogicals(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=mxGetDoubles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=mxGetSingles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=mxGetInt8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=mxGetUint8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=mxGetInt16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=mxGetUint16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=mxGetInt32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=mxGetUint32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=mxGetInt64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=mxGetUint64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *allocUnsignedIntMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(unsigned int)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *convert2DReal2UnsignedIntMat(const mxArray * const val) {
    mxArray *retMat;
    unsigned int *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocUnsignedIntMatInMatlab(M,N);
    if(sizeof(unsigned int)==4) {//32 bit
        retData=(unsigned int*)mxGetUint32s(retMat);
    } else {//64 bit
        retData=(unsigned int*)mxGetUint64s(retMat);
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *data=mxGetChars(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *data=mxGetLogicals(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=mxGetDoubles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=mxGetSingles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=mxGetInt8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=mxGetUint8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=mxGetInt16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=mxGetUint16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=mxGetInt32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=mxGetUint32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=mxGetInt64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=mxGetUint64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *allocSignedSizeMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(ptrdiff_t)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *convert2DReal2SignedSizeMat(const mxArray * const val) {
    mxArray *retMat;
    ptrdiff_t *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocSignedSizeMatInMatlab(M,N);
    if(sizeof(ptrdiff_t)==4) {//32 bit
        retData=(ptrdiff_t*)mxGetInt32s(retMat);
    } else {//64 bit
        retData=(ptrdiff_t*)mxGetInt64s(retMat);
    }

    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *data=mxGetChars(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *data=mxGetLogicals(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=mxGetDoubles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=mxGetSingles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=mxGetInt8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=mxGetUint8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=mxGetInt16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=mxGetUint16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=mxGetInt32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=mxGetUint32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=mxGetInt64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=mxGetUint64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *allocUnsignedSizeMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(size_t)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *convert2DReal2UnsignedSizeMat(const mxArray * const val) {
    mxArray *retMat;
    size_t *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocUnsignedSizeMatInMatlab(M,N);
    if(sizeof(size_t)==4) {//32 bit
        retData=(size_t*)mxGetUint32s(retMat);
    } else {//64 bit
        retData=(size_t*)mxGetUint64s(retMat);
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *data=mxGetChars(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *data=mxGetLogicals(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=mxGetDoubles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=mxGetSingles(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=mxGetInt8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=mxGetUint8s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=mxGetInt16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=mxGetUint16s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=mxGetInt32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=mxGetUint32s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=mxGetInt64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=mxGetUint64s(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

int getIntFromMatlab(const mxArray * const val) {
    int retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(int)*mxGetChars(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(int)*mxGetLogicals(val);
            break;
        case mxDOUBLE_CLASS:
            retVal=(int)*mxGetDoubles(val);
            break;
        case mxSINGLE_CLASS:
            retVal=(int)*mxGetSingles(val);
            break;
        case mxINT8_CLASS:
            retVal=(int)*mxGetInt8s(val);
            break;
        case mxUINT8_CLASS:
            retVal=(int)*mxGetUint8s(val);
            break;
        case mxINT16_CLASS:
            retVal=(int)*mxGetInt16s(val);
            break;
        case mxUINT16_CLASS:
            retVal=(int)*mxGetUint16s(val);
            break;
        case mxINT32_CLASS:
            retVal=(int)*mxGetInt32s(val);
            break;
        case mxUINT32_CLASS:
            retVal=(int)*mxGetUint32s(val);
            break;
        case mxINT64_CLASS:
            retVal=(int)*mxGetInt64s(val);
            break;
        case mxUINT64_CLASS:
            retVal=(int)*mxGetUint64s(val);
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

size_t getSizeTFromMatlab(const mxArray * const val) {
    size_t retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(size_t)*mxGetChars(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(size_t)*mxGetLogicals(val);
            break;
        case mxDOUBLE_CLASS:
        {
            double temp=*mxGetDoubles(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxSINGLE_CLASS:
        {
            float temp=*mxGetSingles(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T temp=*mxGetInt8s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT8_CLASS:
            retVal=(size_t)*mxGetUint8s(val);
            break;
        case mxINT16_CLASS:
        {
            int16_T temp=*mxGetInt16s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT16_CLASS:
            retVal=(size_t)*mxGetUint16s(val);
            break;
        case mxINT32_CLASS:
        {
            int32_T temp=*mxGetInt32s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT32_CLASS:
            retVal=(size_t)*mxGetUint32s(val);
            break;
        case mxINT64_CLASS:
        {
            int64_T temp=*mxGetInt64s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT64_CLASS:
            retVal=(size_t)*mxGetUint64s(val);
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

ptrdiff_t getPtrDiffTFromMatlab(const mxArray * const val) {
    ptrdiff_t retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(ptrdiff_t)*mxGetChars(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(ptrdiff_t)*mxGetLogicals(val);
            break;
        case mxDOUBLE_CLASS:
        {
            double temp=*mxGetDoubles(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxSINGLE_CLASS:
        {
            float temp=*mxGetSingles(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T temp=*mxGetInt8s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxUINT8_CLASS:
            retVal=(ptrdiff_t)*mxGetUint8s(val);
            break;
        case mxINT16_CLASS:
        {
            int16_T temp=*mxGetInt16s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxUINT16_CLASS:
            retVal=(ptrdiff_t)*mxGetUint16s(val);
            break;
        case mxINT32_CLASS:
        {
            int32_T temp=*mxGetInt32s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxUINT32_CLASS:
            retVal=(ptrdiff_t)*mxGetUint32s(val);
            break;
        case mxINT64_CLASS:
        {
            int64_T temp=*mxGetInt64s(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(ptrdiff_t)temp;
            break;
        }
        case mxUINT64_CLASS:
            retVal=(ptrdiff_t)*mxGetUint64s(val);
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

size_t *copySizeTArrayFromMatlab(const mxArray * const val, size_t *arrayLen) {
//The length of the returned array is placed into the arrayLen variable.
    size_t M,N,i;
    size_t numEl=0;
    size_t *retVal=NULL;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfDimensions(val)>2) {
        mexErrMsgTxt("A parameter that should be an array has extra dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    
    if(M==1&&N>=1) {
        numEl=N;
    } else if(M>1&&N==1){
        numEl=M;
    } else {
        mexErrMsgTxt("A parameter that should be an array has an incorrect dimensionality.");
    }
    
    *arrayLen=numEl;
    retVal=(size_t*)mxMalloc(numEl*sizeof(size_t));    
            
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *mexData=mxGetChars(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *mexData=mxGetLogicals(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *mexData=mxGetDoubles(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *mexData=mxGetSingles(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *mexData=mxGetInt8s(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *mexData=mxGetUint8s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *mexData=mxGetInt16s(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *mexData=mxGetUint16s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *mexData=mxGetInt32s(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *mexData=mxGetUint32s(val);

            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *mexData=mxGetInt64s(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    mxFree(retVal);
                    
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *mexData=mxGetUint64s(val);

            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mxFree(retVal);
            
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retVal;
}

ptrdiff_t *copyPtrDiffTArrayFromMatlab(const mxArray * const val, size_t *arrayLen) {
//The length of the returned array is placed into the arrayLen variable.
    size_t M,N,i;
    size_t numEl=0;
    ptrdiff_t *retVal=NULL;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfDimensions(val)>2) {
        mexErrMsgTxt("A parameter that should be an array has extra dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    
    if(M==1&&N>=1) {
        numEl=N;
    } else if(M>1&&N==1){
        numEl=M;
    } else {
        mexErrMsgTxt("A parameter that should be an array has an incorrect dimensionality.");
    }
    
    *arrayLen=numEl;
    retVal=(ptrdiff_t*)mxMalloc(numEl*sizeof(ptrdiff_t));    
            
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *mexData=mxGetChars(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *mexData=mxGetLogicals(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *mexData=mxGetDoubles(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *mexData=mxGetSingles(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *mexData=mxGetInt8s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *mexData=mxGetUint8s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *mexData=mxGetInt16s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *mexData=mxGetUint16s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *mexData=mxGetInt32s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *mexData=mxGetUint32s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *mexData=mxGetInt64s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *mexData=mxGetUint64s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(ptrdiff_t)mexData[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mxFree(retVal);
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retVal;
}

bool *copyBoolArrayFromMatlab(const mxArray * const val, size_t *arrayLen) {
//The length of the returned array is placed into the arrayLen variable.
    size_t M,N,i;
    size_t numEl=0;
    bool *retVal=NULL;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfDimensions(val)>2) {
        mexErrMsgTxt("A parameter that should be an array has extra dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    
    if(M==1&&N>=1) {
        numEl=N;
    } else if(M>1&&N==1){
        numEl=M;
    } else {
        mexErrMsgTxt("A parameter that should be an array has an incorrect dimensionality.");
    }
    
    *arrayLen=numEl;
    retVal=(bool*)mxMalloc(numEl*sizeof(bool));  
            
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            mxChar *mexData=mxGetChars(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical *mexData=mxGetLogicals(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *mexData=mxGetDoubles(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *mexData=mxGetSingles(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *mexData=mxGetInt8s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *mexData=mxGetUint8s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *mexData=mxGetInt16s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *mexData=mxGetUint16s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *mexData=mxGetInt32s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *mexData=mxGetUint32s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *mexData=mxGetInt64s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *mexData=mxGetUint64s(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(bool)mexData[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mxFree(retVal);
            
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retVal;
}

double getDoubleFromMatlab(const mxArray * const val) {
    double retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(double)*mxGetChars(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(double)*mxGetLogicals(val);
            break;
        case mxDOUBLE_CLASS:
            retVal=(double)*mxGetDoubles(val);
            break;
        case mxSINGLE_CLASS:
            retVal=(double)*mxGetSingles(val);
            break;
        case mxINT8_CLASS:
            retVal=(double)*mxGetInt8s(val);
            break;
        case mxUINT8_CLASS:
            retVal=(double)*mxGetUint8s(val);
            break;
        case mxINT16_CLASS:
            retVal=(double)*mxGetInt16s(val);
            break;
        case mxUINT16_CLASS:
            retVal=(double)*mxGetUint16s(val);
            break;
        case mxINT32_CLASS:
            retVal=(double)*mxGetInt32s(val);
            break;
        case mxUINT32_CLASS:
            retVal=(double)*mxGetUint32s(val);
            break;
        case mxINT64_CLASS:
            retVal=(double)*mxGetInt64s(val);
            break;
        case mxUINT64_CLASS:
            retVal=(double)*mxGetUint64s(val);
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

bool getBoolFromMatlab(const mxArray * const val) {
    bool retVal=false;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=*mxGetChars(val)!=0.0;
            break;
        case mxLOGICAL_CLASS:
            retVal=*mxGetLogicals(val)!=0;
            break;
        case mxDOUBLE_CLASS:
            retVal=*mxGetDoubles(val)!=0.0;
            break;
        case mxSINGLE_CLASS:
            retVal=*mxGetSingles(val)!=0;
            break;
        case mxINT8_CLASS:
            retVal=*mxGetInt8s(val)!=0;
            break;
        case mxUINT8_CLASS:
            retVal=*mxGetUint8s(val)!=0;
            break;
        case mxINT16_CLASS:
            retVal=*mxGetInt16s(val)!=0;
            break;
        case mxUINT16_CLASS:
            retVal=*mxGetUint16s(val)!=0;
            break;
        case mxINT32_CLASS:
            retVal=*mxGetInt32s(val)!=0;
            break;
        case mxUINT32_CLASS:
            retVal=*mxGetUint32s(val)!=0;
            break;
        case mxINT64_CLASS:
            retVal=*mxGetInt64s(val)!=0;
            break;
        case mxUINT64_CLASS:
            retVal=*mxGetUint64s(val)!=0;
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

int32_t getInt32FromMatlab(const mxArray * const val) {
    int32_t retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(int32_t)*mxGetChars(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(int32_t)*mxGetLogicals(val);
            break;
        case mxDOUBLE_CLASS:
            retVal=(int32_t)*mxGetDoubles(val);
            break;
        case mxSINGLE_CLASS:
            retVal=(int32_t)*mxGetSingles(val);
            break;
        case mxINT8_CLASS:
            retVal=(int32_t)*mxGetInt8s(val);
            break;
        case mxUINT8_CLASS:
            retVal=(int32_t)*mxGetUint8s(val);
            break;
        case mxINT16_CLASS:
            retVal=(int32_t)*mxGetInt16s(val);
            break;
        case mxUINT16_CLASS:
            retVal=(int32_t)*mxGetUint16s(val);
            break;
        case mxINT32_CLASS:
            retVal=(int32_t)*mxGetInt32s(val);
            break;
        case mxUINT32_CLASS:
            retVal=(int32_t)*mxGetUint32s(val);
            break;
        case mxINT64_CLASS:
            retVal=(int32_t)*mxGetInt64s(val);
            break;
        case mxUINT64_CLASS:
            retVal=(int32_t)*mxGetUint64s(val);
            break;
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal; 
}

mxArray *doubleMat2Matlab(const double * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    double *dataPtr;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    memcpy(dataPtr,arr,numRow*numCol*sizeof(double));

    return retMat;
}

mxArray *floatMat2MatlabDoubles(const float * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    double *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=(double)arr[i];
    }

    return retMat;
}

mxArray *intMat2MatlabDoubles(const int * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    double *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=(double)arr[i];
    }

    return retMat;
}

mxArray *sizeTMat2MatlabDoubles(const size_t * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    double *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=(double)arr[i];
    }

    return retMat;
}

mxArray *ptrDiffTMat2MatlabDoubles(const ptrdiff_t * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    double *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=(double)arr[i];
    }

    return retMat;
}

mxArray *boolMat2Matlab(const bool * const arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    mxLogical *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateLogicalMatrix(numRow,numCol);
    dataPtr=mxGetLogicals(retMat);
    
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=(mxLogical)arr[i];
    }

    return retMat;
}

mxArray *allocUnsignedCharMatInMatlab(const size_t numRow, const size_t numCol) {
    return mxCreateNumericMatrix(numRow,numCol,mxUINT8_CLASS,mxREAL);
}

mxArray *unsignedSizeMat2Matlab(const size_t * const arr, const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    size_t *dataPtr;
    
    retMat=allocUnsignedSizeMatInMatlab(numRow, numCol);
    
    if(sizeof(size_t)==4) {//32 bit
        dataPtr=(size_t*)mxGetUint32s(retMat);
    } else {//64 bit
        dataPtr=(size_t*)mxGetUint64s(retMat);
    }

    memcpy(dataPtr,arr,numRow*numCol*sizeof(size_t));

    return retMat;
}

mxArray *unsignedCharMat2Matlab(const unsigned char * const arr, const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    unsigned char *dataPtr;
    
    retMat=allocUnsignedCharMatInMatlab(numRow, numCol);
    
    dataPtr=(unsigned char*)mxGetUint8s(retMat);
    memcpy(dataPtr,arr,numRow*numCol*sizeof(unsigned char));

    return retMat;
}

mxArray *signedSizeMat2Matlab(const ptrdiff_t * const arr, const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    ptrdiff_t *dataPtr;
    
    retMat=allocSignedSizeMatInMatlab(numRow,numCol);
    
    if(sizeof(ptrdiff_t)==4) {//32 bit
        dataPtr=(ptrdiff_t*)mxGetInt32s(retMat);
    } else {//64 bit
        dataPtr=(ptrdiff_t*)mxGetInt64s(retMat);
    }

    memcpy(dataPtr,arr,numRow*numCol*sizeof(ptrdiff_t));

    return retMat;
}

double getScalarMatlabClassConst(const char * const className, const char * const constName) {
/*Whereas in Matlab one can access constants in classes without
 *instantiating the class, in mex, the class must be instantiated first.
 *Class name and const name are null-terminated strings. For example, if one wishes to access
 *Constants.standardTemp, then one can use
 *const char className[]="Constants";
 *const char constName[]="standardTemp";
 **/
    mxArray *constantClass, *valMATLAB;
    double retVal;

    mexCallMATLAB(1,&constantClass,0,NULL,className);//Load the Constants class.
    valMATLAB=mxGetProperty(constantClass, 0,constName);

    if(valMATLAB==NULL) {
        mexErrMsgTxt("An eror occurred while trying to access a class constant.");
    }
    mxDestroyArray(constantClass);
    retVal=getDoubleFromMatlab(valMATLAB);
    mxDestroyArray(valMATLAB);
    return retVal;
}

bool pointerIsAligned(const void *pointer, const size_t byteCount) {
    return (uintptr_t)pointer%byteCount == 0;
}

#ifdef __cplusplus

template<typename T>
mxArray *ptr2Matlab(T thePointer) {
/*Ptr2Matlab Convert a C++ pointer to a Matlab matrix so that it can be
 *           returned to Matlab between calls to mex functions.
 */
    mxArray *retArray=NULL;
    
    switch(sizeof(void*)) {
        case 4:
            retArray=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retArray=mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    *reinterpret_cast<void**>(mxGetData(retArray))=reinterpret_cast<void*>(thePointer);
    return retArray;
}

template<typename T>
T Matlab2Ptr(const mxArray * const matlabPtr) {
/**MATLAB2PTR Convert a 1X1 Matlab matrix in which a pointer has been saved
 *            using ptr2Matlab back into a pointer for use in a mex file.
 */
    T thePtr;
    thePtr=reinterpret_cast<T>(*reinterpret_cast<void**>(mxGetData(matlabPtr)));    
    return thePtr;
}


template<typename T>
mxArray *mat2MatlabDoubles(T *arr,const size_t numRow,const size_t numCol) {
/**MAT2MATLABDOUBLES Convert a matrix of some template format into a matrix
 *                   of doubles to be returned to Matlab by copying each
 *                   element into the Matlab matrix.
 */
    mxArray *retMat;
    double *dataPtr;
    const size_t totalNumEl=numRow*numCol;
    size_t i;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetDoubles(retMat);
    
    for(i=0;i<totalNumEl;i++) {
        dataPtr[i]=static_cast<double>(arr[i]);   
    }

    return retMat;
}

std::complex <double> getComplexDoubleFromMatlab(const mxArray * const val) {

    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    if(mxIsComplex(val)==false) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
 
    switch(mxGetClassID(val)){
        case mxDOUBLE_CLASS:
        {
            mxComplexDouble *theCompNum=mxGetComplexDoubles(val);
            return std::complex<double>(theCompNum->real,theCompNum->imag);
        }
        case mxSINGLE_CLASS:
        {
            mxComplexSingle *theCompNum=mxGetComplexSingles(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxINT8_CLASS:
        {
            mxComplexInt8 *theCompNum=mxGetComplexInt8s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxUINT8_CLASS:
        {
            mxComplexUint8 *theCompNum=mxGetComplexUint8s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxINT16_CLASS:
        {
            mxComplexInt16 *theCompNum=mxGetComplexInt16s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxUINT16_CLASS:
        {
            mxComplexUint16 *theCompNum=mxGetComplexUint16s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxINT32_CLASS:
        {
            mxComplexInt32 *theCompNum=mxGetComplexInt32s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxUINT32_CLASS:
        {
            mxComplexUint32 *theCompNum=mxGetComplexUint32s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxINT64_CLASS:
        {
            mxComplexInt64 *theCompNum=mxGetComplexInt64s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxUINT64_CLASS:
        {
            mxComplexUint64 *theCompNum=mxGetComplexUint64s(val);
            return std::complex<double> (static_cast<double>(theCompNum->real),static_cast<double>(theCompNum->imag));
        }
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("A parameter is of a data type that cannot be used.");
    }
    
    //It should never get here.
    return 0;
}

#endif
#endif

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
*/
