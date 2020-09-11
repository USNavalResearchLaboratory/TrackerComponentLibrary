/**EXACTSIGNOFSUMCPP The main C++ routine for computing the exact sign of
 *                   the sum of a number of finite numbers. The sign
 *                   computation occurs exactly, regardless of whether the
 *                   actual sum value would have overflowed or suffered an
 *                   accumulation of finite precision errors whereby just
 *                   adding them would have given the wrong result.
 *
 *The entire implementation is given in the header file so that it can be a
 *template class. This is defined as a template class so that it can be
 *used with float as well as double and long double values. It cannot be
 *used with integers.
 *
 *Only the function int exactSignOfSumCPP is meant to be called externally.
 *
 *The algorithm is taken from [2] where code is provided in an appendix.
 *The code with minor changes and corrections is also available from Jon
 *Rokne's web site at
 *http://pages.cpsc.ucalgary.ca/~rokne/convex/sgnsum.cc
 *The implementation here uses the corrections and uses the sort algorithm
 *in the C++ standard template library rather than the sort algorithm
 *provided by Rokne.
 *
 *REFERENCES:
 *[1] H. Ratschek and J. Rokne, "Exact computation of the sign of a finite
 *   sum," Applied Mathematics and Computation, vol. 99, no. 2-3, pp.
 *   99-127, 15 Mar. 1999.
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef EXACTSIGNOFSUMCPP
#define EXACTSIGNOFSUMCPP

//Needed for sort
#include <algorithm>
//Needed for greater
#include <functional>
//If C++11 with type-generic versions of frexp and ldexp are supported.
#if __cplusplus<=199711L
//If no type generic version of frexp is available, then that means that
//we just have to use cmath and deal with an increased possibility of
//overflows
#include <cmath>
#else
//Needed for type-generic frexp and ldexp
#include <ctgmath>
#endif

//Two helper functions for keeping things in order.
template<typename T>
void BuildHeapFromTop(const size_t n, T *ra)
{
 size_t i=1,m;
 T top=ra[1];

 while (2*i<=n)
 {
  m= 2*i;
  if (ra[m]<ra[m+1]) if (m<n) m++;
  if (top<ra[m]) { ra[i]=ra[m]; i=m; }
  else break; 
 }
 ra[i]=top;
}

template<typename T>
void BuildHeapFromBelow(const size_t n, T *ra)
{
 size_t i=n,m;
 T last=ra[n];
 
 while (i/2>0)
 {
  m= i/2;
  if ( ra[m]<last ) { ra[i]=ra[m]; i=m;}
  else break;
 }
 ra[i]= last;
}

template<typename T>
int exactSignOfSumCPP(const T *S, const size_t nS) {
    int sg;//The return value
    size_t n, m, i;
    T *a, *b;
    //Initialization of the lists a and b. The extra parentheses zeros the
    //array of doubles/ floats so that the call to new behaves like calloc.
    a=new T[nS+3]();
    b=new T[nS+3]();
    //Splitting of S into positive and negative summands.
    n=0;
    m=0;
    for(i=0;i<nS;i++) {
        if(S[i]>0){
            a[++m]=S[i];
        } else if(S[i]<0) {
            b[++n]=-S[i];
        }
    }
    
    //Building the lists a and b into heaps descending order
    if(m>1) {
        std::sort(a+1,a+m+1,std::greater<T>());
    }
    if(n>1) {
        std::sort(b+1,b+n+1,std::greater<T>());
    }

    {
        int E, F;
        T as, ass;
        T bs, bss;
        T uu, u, v;
        // Main loop (the proper ESSA algorithm)
        
        while(1) {
            //Step 1: (Termination Criteria)
            if (n==0&&m==0){
                sg=0;
                break;
            } else if (n==0||(a[1]>n*b[1]&&m>0)){
                sg=1;
                break;
            } else if (m==0||(b[1]>m*a[1]&&n>0)){
                sg=-1;
                break;
            }
            
            //Step 2: (Auxiliary variables)
            as=ass=bs=bss=0;
            
            //Step 3: (Compare and process the leading summands of the
            //lists
            
            //Computation of the exponent E from a[1] and F from b[1] with
            //respect to base 2. a[1]-2^E (b[1]-2^F) zeros the first bit of
            //a[1] (b[1])
            frexp(a[1],&E);
            E--;
            frexp(b[1],&F);
            F--;

            //E contains the exponent of a[1] and F contains the exponent
            //of b[1] in base 2.
            if(E==F)//Step 3, case (i):
            {
                if(a[1]>=b[1]){
                    as=a[1]-b[1];
                } else{
                    bs=b[1]-a[1];
                }
            }
            else if(E>F)//Step 3, case (ii):
            {
                uu=ldexp(1.0,F);
                if(b[1]==uu) {
                    u=uu;
                } else {
                    u=2*uu;
                }

                as=a[1]-u;
                ass=u-b[1];
            }
            else if(F>E)//Step 3, case (iii):
            {
                uu=ldexp(1.0,E);
                if(a[1]==uu) {
                    v=uu;
                } else {
                    v=2*uu;
                }
                
                bs= b[1]-v;
                bss=v-a[1];
            }
            
            // Step 4: Rearrangement of the lists, keeping S constant.
            if(as==0) {
                if(ass==0) {
                    a[1]=a[m];
                    m--;
                } else {
                    a[1]=ass;
                }
                BuildHeapFromTop(m,a);
            } else {
                a[1]=as;
                BuildHeapFromTop(m,a);
                if(ass!=0) {
                    a[++m]=ass;   
                    BuildHeapFromBelow(m,a);
                }
            }
            
            if(bs==0) {
                if(bss==0) {
                    b[1]=b[n];
                    n--;
                } else {
                    b[1]=bss;
                }
                BuildHeapFromTop(n,b);
            } else {
                b[1]=bs;
                BuildHeapFromTop(n,b);
                if(bss!=0) {
                    b[++n]=bss;
                    BuildHeapFromBelow(n,b);
                }
            }
        }
        
        delete a;
        delete b;
        
        return sg;
    }
}

#endif

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
