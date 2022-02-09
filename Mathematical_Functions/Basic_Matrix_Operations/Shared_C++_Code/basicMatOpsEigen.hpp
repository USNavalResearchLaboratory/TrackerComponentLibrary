/**BASICMATOPSEIGEN C++ functions of basic operations with matrices, having
 *             Eigen as a dependency.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef BASICMATOPSEIGENCPP
#define BASICMATOPSEIGENCPP

//Has isfinite, NAN, isnan and nextafter
#include <cmath>
//Defines the size_t
#include <cstddef>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/Eigenvalues"
#include <complex>
//For infinity
#include <limits>
//For max
#include <algorithm>

/**ANYNONFINITEINMATRIX Given an Eigen matrix type, this returns true is
 *          any non-finite values are therein. Otherwise, it returns false.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
bool anyNonFiniteInMatrix(const T &theMat) {
    const size_t numRow=theMat.rows();
    const size_t numCol=theMat.cols();
    const size_t numEls=numRow*numCol;
    
    for(size_t k=0;k<numEls;k++) {
        if(!std::isfinite(std::abs(theMat(k)))) {
            return true;
        }
    }
    return false;
}

/**MATRIXOFVAL Return an Eigen matrix of the specified dimensions where all
 *             elements are the same value, val. It is assumed that the
 *             template type TEig is something like Eigen::MatrixXd or
 *             Eigen::MatrixXcd which would mean that TVal should be
 *             something like double or std::complex<double>.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class TEig, class TVal>
TEig MatrixOfVal(const TVal val,const size_t numRows, const size_t numCols) {
    TEig theMat(numRows,numCols);
    const size_t numEls=numRows*numCols;
    
    for(size_t k=0;k<numEls;k++) {
       theMat(k)=val;
    }
      
    return theMat;
}

/**MATRIXRANK Determine the rank of a matrix using a chosen criterion. This
 *       offers more flexibility than the rank option in Eigen. This
 *       function works with input matrices of type Eigen::MatrixXd and
 *       type Eigen::MatrixXcd. 
 *
 *INPUTS: X An MXN matrix whose rank is desired.
 *  algorithm This specifies the criterion for determining the rank.
 *          0 Use the Matlab equivalent of eps(norm(X,1)) as the bound.
 *            That means, this is the difference between the smallest
 *            double floating point number larger than norm(X,1) and
 *            norm(X,1) (the 1-norm fo the matrix).
 *          1 Use the Matlab equivalent of max(size(X))*eps(max(s)), where
 *            s is the vector of singular values. eps(max(s)) is the
 *            floating point difference between the next largest number 
 *            after max(s) and max(s).
 *          2 Use std::numeric_limits<double>::epsilon()*norm(A,1).
 *          3 Use max(size(A))*std::numeric_limits<double>::epsilon()*max(s).
 *    U,s,V These are returned as the left (U) and right (V) singular
 *          vectors (not computed thin) and a vector of the singular
 *          values. 
 *
 *OUTPUTS: The return value is the rank. The outputs of the SVD performed
 *         in computing the rank are put in U, s, and V for reuse
 *         elsewhere.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. Thus, the function might be called as
 *matrixRank<Eigen::MatrixXd>(X, algorithm,U, s, V)
 *
 *The notion of the matrix rank relating to a thresholding of the singular
 *values is discussed in Chapter 5.4.1 of [1].
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
size_t matrixRank(const T &X, const int algorithm,T &U, Eigen::MatrixXd &s, T &V) {
    const size_t M=X.rows();
    const size_t N=X.cols();
    const size_t sLen=std::min(M,N);
    
    //Check for non-finite values. If any are found, return a bunch of NaNs
    //in U, s and V and a zero rank.
    if(anyNonFiniteInMatrix(X)) {
        U=MatrixOfVal<T,double>(NAN,M,sLen);
        V=MatrixOfVal<T,double>(NAN,N,sLen);
        s=MatrixOfVal<Eigen::MatrixXd,double>(NAN,sLen,1);

        return 0;
    }
    
    Eigen::JacobiSVD<T> svdX(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //The singular values are already sorted in decreasing order.
    s=svdX.singularValues();
    U=svdX.matrixU();
    V=svdX.matrixV();

    const double epsVal=std::numeric_limits<double>::epsilon();
    double limVal;
    switch(algorithm) {
        case 0:
        {
            //Get the 1-matrix norm of X. This is the maximum absolute
            //column sum of the matrix.
            double normX=X.col(0).cwiseAbs().sum();
            
            for(size_t i=1;i<N;i++) {
                const double curAbsColSum=X.col(i).cwiseAbs().sum();
                normX=std::max(curAbsColSum,normX);
            }

            limVal=nextafter(normX,std::numeric_limits<double>::infinity())-normX;
            break;
        }
        case 1:
        {
            limVal=static_cast<double>(std::max(N,M))*(nextafter(s(0),std::numeric_limits<double>::infinity())-s(0));
            break;
        }
        case 2:
        {
            //Get the 1-matrix norm of X. This is the maximum absolute
            //column sum of the matrix.
            double normX=X.col(0).cwiseAbs().sum();
            for(size_t i=1;i<N;i++) {
                const double curAbsColSum=X.col(i).cwiseAbs().sum();
                normX=std::max(curAbsColSum,normX);
            }
            limVal=epsVal*normX;
            break;
        }
        default:
        {
            limVal=static_cast<double>(std::max(N,M))*s(0)*epsVal;
            break;
        }
    };
    
    size_t rankVal=0;
    for(size_t k=0;k<sLen;k++) {
        if(s(k)>limVal) {
            rankVal++;
        } else {
            break;
        }
    }
    
    return rankVal;
}

/**PSEUDOINVERSE Compute the matrix pseudoinverse using a singular value
 *       decomposition. This is essentially equivalent to the pinv function
 *       in Matlab, except one can choose the algorithm to use for
 *       determining the rank of the matrix. This function works with input
 *       matrices of type Eigen::MatrixXd and type Eigen::MatrixXcd. 
 *
 *INPUTS: X An MXN matrix whose pseudoinverse is desired.
 *  algorithm This specifies the criterion for determining the rank, which
 *          is a necessary step in computing the pseudoinverse.
 *          0 Use the Matlab equivalent of eps(norm(X,1)) as the bound.
 *            That means, this is the difference between the smallest
 *            double floating point number larger than norm(X,1) and
 *            norm(X,1) (the 1-norm fo the matrix).
 *          1 Use the Matlab equivalent of max(size(X))*eps(max(s)), where
 *            s is the vector of singular values. eps(max(s)) is the
 *            floating point difference between the next largest number 
 *            after max(s) and max(s).
 *          2 Use std::numeric_limits<double>::epsilon()*norm(A,1).
 *          3 Use max(size(A))*std::numeric_limits<double>::epsilon()*max(s).
 *
 *OUTPUTS: The return value is the NXM pseudoinverse.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. Thus, the function might be called as
 *pseudoInverse<Eigen::MatrixXd>(X, algorithm,U, s, V)
 *
 *The notion of the matrix rank relating to a thresholding of the singular
 *values is discussed in Chapter 5.4.1 of [1].
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
T pseudoInverse(const T X,const int algorithm) {
    T U;
    Eigen::MatrixXd s;
    T V;
    
    const size_t rankVal=matrixRank(X, algorithm, U, s, V);
    
    //The singular values in S are in descending order.
    const size_t N=s.rows();
    for(size_t k=0;k<rankVal;k++) {
        s(k)=1/s(k);
    }
    for(size_t k=rankVal;k<N;k++) {
        s(k)=0.0;
    }
    
    T retVal=V*s.asDiagonal()*U.transpose();
    return retVal;
}

/**FORWARDSUBSTITUTIONBYROW Given a lower-triangular matrix, use forward
 *      substitution to solve the system L*x=b for x. This function uses
 *      the row-oriented algorithm, Algorithm 3.1.1 in [1].
 *
 *INPUTS: L An NXN lower-triangular matrix.
 *        b An NX1 vector in which the result is also placed.
 *
 *OUTPUTS: The output is put in b.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. For example, one might call the function as
 *forwardSubstitutionByRow<Eigen::MatrixXd>(L,b,algorithm)
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void forwardSubstitutionByRow(const T &L, T &b) {
    const size_t n=L.rows();
    b(0)=b(0)/L(0,0);
    for(size_t k=1;k<n;k++) {
        double sumVal=0;
        for(size_t i=0;i<k;i++) {
            sumVal+=L(k,i)*b(i);
        }        
        b(k)=(b(k)-sumVal)/L(k,k);
    }
}

/**FORWARDSUBSTITUTIONBYCOL Given a lower-triangular matrix, use forward
 *      substitution to solve the system L*x=b for x. This function uses
 *      the column-oriented algorithm, Algorithm 3.1.3 in [1].
 *
 *INPUTS: L An NXN lower-triangular matrix.
 *        b An NX1 vector in which the result is also placed.
 *
 *OUTPUTS: The output is put in b.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. For example, one might call the function as
 *forwardSubstitutionByCol<Eigen::MatrixXd>(L,b,algorithm)
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void forwardSubstitutionByCol(const T &L, T &b) {
    const size_t n=L.rows();
    //Algorithm 3.1.3 in [1], column-oriented.
    for(size_t k=0;k<n-1;k++) {
        b(k)=b(k)/L(k,k);
        for(size_t i=k+1;i<n;i++) {
            b(i)-=L(i,k)*b(k);
        }
    }          
    b(n-1)=b(n-1)/L(n-1,n-1);
}

/**BACKSUBSTITUTIONBYROW Given an upper-triangular matrix, use back
 *      substitution to solve the system U*x=b for x.  This function uses
 *      the row-oriented algorithm, Algorithm 3.1.2 in in [1].
 *
 *INPUTS: U An NXN upper-triangular matrix.
 *        b An NX1 vector in which the result is also placed.
 *
 *OUTPUTS: The output is put in b.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. For example, one might call the function as
 *backSubstitutionByRow<Eigen::MatrixXd>(L,b,algorithm)
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void backSubstitutionByRow(const T &U, T &b) {
    const size_t n=U.rows();
    //Algorithm 3.1.2 in [1], row-oriented.
    b(n-1)=b(n-1)/U(n-1,n-1);
    size_t k=n-1;
    do {
        k--;
        double sumVal=0;
        for(size_t i=k+1;i<n;i++) {
            sumVal+=U(k,i)*b(i);
        }
        b(k)=(b(k)-sumVal)/U(k,k);
    } while(k>0);
}

/**BACKSUBSTITUTIONBYCOL Given an upper-triangular matrix, use back
 *      substitution to solve the system U*x=b for x.  This function uses
 *      the columns-oriented algorithm, Algorithm 3.1.4 in in [1].
 *
 *INPUTS: U An NXN upper-triangular matrix.
 *        b An NX1 vector in which the result is also placed.
 *
 *OUTPUTS: The output is put in b.
 *
 *T is an Eigen class. The class often needs to be explicitly given when
 *using this function. For example, one might call the function as
 *backSubstitutionByCol<Eigen::MatrixXd>(L,b,algorithm)
 *
 *REFERENCES:
 *[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
 *   Baltimore: Johns Hopkins University Press, 2013.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void backSubstitutionByCol(const T &U, T &b) {
    const size_t n=U.rows();
    
    for(size_t k=n-1;k>=1;k--) {
        b(k)=b(k)/U(k,k);
        for(size_t i=0;i<k;i++) {
            b(i)-=U(i,k)*b(k);
        }
    }
    b(0)=b(0)/U(0);
}

/*TWOMATDIAGSVD Return an NXN matrix that diagonalizes two real or complex
 *              matrices. This is found using an SVD-based algorithm. Both
 *              inputs must be either of type Eigen::MatrixXd or of type
 *              Eigen::MatrixXcd. 
 *
 *INPUTS: C1, C2 The two NXN complex matrices to diagonalize. These are
 *               stored by column.
 *
 *OUTPUTS: An NXN matrix of the same type as the input is returned.
 *
 *Note that in JacobiSVD in some versions of Eigen, the return arrays are
 *not fully allocated if NaNs or Inf values are present. Thus, to avoid
 *reading past the end of an array, if NaNs or Infs arise, this function
 *returns all NaNs.
 *
 *Note, normally with the function below, one needs to explicitly list the
 *class. For example, twoMatDiagSVD<MatrixXd>(C1,C2) or
 *twoMatDiagSVD<MatrixXcd>(C1,C2).
 *
 *The algorithm of [1] is used.
 *
 *REFERENCES:
 *[1] J. Nygårds, V. Deleskog, and G. Hendeby, "Safe fusion compared to
 *    established distributed fusion methods," in IEEE International
 *    Conference on Multisensor Fusion and Integration for Intelligent
 *    Systems, Baden-Baden, Germany, 19-21 Sep. 2016, pp. 265-271.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template <class T>
T twoMatDiagSVD(const T &C1,const T &C2) {
    const size_t N=C1.rows();
    const size_t N2=N*N;

    //Check that we are not passing JacobiSVD any NaNs or Inf values. If we
    //are, return a matrix of NaNs without calling JacobiSVD. This is done
    //explicitly, because newer versions of Eigen set an info parameter
    //marking the it as invalid, but do not actually compute U and can
    //have invalid memory access issues if one tries to use the result. We
    //must do this for all calls to JacobiSVD. We do not directly call
    //isfinite, because it won't work if the type here is complex. Thus, we
    //first take the absolute value of the input.
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(std::abs(C1(k1)))) {
            T retMat(N,N);
            
            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }
            return retMat;
        }
    }

    const Eigen::JacobiSVD<T> C1SVD(C1, Eigen::ComputeFullU);
    T d=C1SVD.singularValues();
    const T U1=C1SVD.matrixU();
    
    d=d.cwiseSqrt().cwiseInverse().eval();
    const T temp=U1*d.asDiagonal();
    T inputMat=temp.adjoint()*C2*temp;
    
    //Again, check that the input to JacobiSVD is valid.
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(std::abs(inputMat(k1)))) {
            T retMat(N,N);
            
            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }
            return retMat;
        }
    }

    const Eigen::JacobiSVD<T> SVD2(inputMat, Eigen::ComputeFullU);
    const T U2=SVD2.matrixU();
    
    //Equation 8a in [1].
    T retVal=U2.adjoint()*d.asDiagonal()*U1.adjoint();
    return retVal;
}

/*TWOMATDIAGEIG Return an NXN matrix that diagonalizes two real or complex
 *              matrices. This is found using an SVD-based algorithm. Both
 *              inputs must be either of type Eigen::MatrixXd or of type
 *              Eigen::MatrixXcd. 
 *
 *INPUTS: C1, C2 The two NXN complex matrices to diagonalize. These are
 *               stored by column.
 *
 *OUTPUTS: An NXN matrix of the same type as the input is returned.
 *
 *Note, normally with the function below, one needs to explicitly list the
 *class. For example, twoMatDiagEig<MatrixXd>(C1,C2) or
 *twoMatDiagEig<MatrixXcd>(C1,C2).
 *
 *The algorithm of [1] is used.
 *
 *REFERENCES:
 *[1] M. Reinhardt, B. Noack, and U. D. Hanebeck, "Closed-form optimization
 *   of covariance intersection for low-dimensional matrices," in 
 *   Proceedings of the 15th International Conference on Information
 *   Fusion, Singapore, 9-12 Jun. 2012, pp. 1891-1896.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template <class T>
T twoMatDiagEig(const T &C1,const T &C2) {};
template <>
Eigen::MatrixXcd twoMatDiagEig(const Eigen::MatrixXcd &C1,const Eigen::MatrixXcd &C2) {
    //Use the eigenvalue-based algorithm.
    const size_t N=C1.rows();
    const size_t N2=N*N;
    
    //Check that we are not passing ComplexEigenSolver any NaNs or Inf
    //values. If we are, return a matrix of NaNs without calling
    //ComplexEigenSolver. This is done explicitly, because newer versions
    //of Eigen set an info parameter marking the it as invalid, but do not
    //actually compute the desired results and can have invalid memory
    //access issues if one tries to use the result. We must do this for all
    //calls to ComplexEigenSolver. isfinite must be given the absolute
    //value of the quantity being tested, because it does not 
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(std::abs(C1(k1)))) {
            Eigen::MatrixXcd retMat(N,N);

            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }

            return retMat;
        }
    }

    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigenC1(C1);
    const Eigen::MatrixXcd V1=eigenC1.eigenvectors();
    Eigen::MatrixXcd d=eigenC1.eigenvalues();
 
    d=d.cwiseSqrt().eval();
    Eigen::MatrixXcd T1=V1*d.asDiagonal();
    T1=T1.inverse().eval();
    const Eigen::MatrixXcd inputMat=T1*C2*T1.adjoint();
    
    //Testing the input before using ComplexEigenSolver.
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(std::abs(C1(k1)))) {
            Eigen::MatrixXcd retMat(N,N);

            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }

            return retMat;
        }
    }
    
    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigenC(inputMat);
    const Eigen::MatrixXcd V2p=eigenC.eigenvectors();
    
    d=d.cwiseInverse().eval();
    //The transformation matrix in Equation 2 of [2].
    Eigen::MatrixXcd retMat=V2p.adjoint()*d.asDiagonal()*V1.adjoint();
    return retMat;
}

template <>
Eigen::MatrixXd twoMatDiagEig(const Eigen::MatrixXd &C1,const Eigen::MatrixXd &C2) {
    //Use the eigenvalue-based algorithm.
    const size_t N=C1.rows();
    const size_t N2=N*N;
    
    //Check that we are not passing EigenSolver any NaNs or Inf values. If
    //we are, return a matrix of NaNs without calling EigenSolver. This is
    //done explicitly, because newer versions of Eigen set an info
    //parameter marking the it as invalid, but do not actually compute the
    //desired results and can have invalid memory access issues if one
    //tries to use the result. We must do this for all calls to
    //EigenSolver.
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(C1(k1))) {
            Eigen::MatrixXd retMat(N,N);

            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }

            return retMat;
        }
    }

    const Eigen::EigenSolver<Eigen::MatrixXd> eigenC1(C1);
    const Eigen::MatrixXd V1=eigenC1.eigenvectors().real();
    Eigen::MatrixXd d=eigenC1.eigenvalues().real();

    d=d.cwiseSqrt().eval();
    Eigen::MatrixXd T1=V1*d.asDiagonal();
    T1=T1.inverse().eval();
    const Eigen::MatrixXd inputMat=T1*C2*T1.adjoint();
    
    for(size_t k1=0;k1<N2;k1++) {
        if(!std::isfinite(inputMat(k1))) {
            Eigen::MatrixXd retMat(N,N);

            for(size_t i1=0;i1<N2;i1++) {
                retMat(i1)=NAN;
            }

            return retMat;
        }
    }
    
    const Eigen::EigenSolver<Eigen::MatrixXd> eigenC(inputMat);
    const Eigen::MatrixXd V2p=eigenC.eigenvectors().real();
    
    d=d.cwiseInverse().eval();
    //The transformation matrix in Equation 2 of [2].
    Eigen::MatrixXd retMat=V2p.adjoint()*d.asDiagonal()*V1.adjoint();
    return retMat;
}

/**COMPUTEFANDDCOSTS This function is used internally to the
 *          JointMatDiagZieheAlg and JointMatDiagFrobVollAlg functions for
 *          computing the costs associated with the matrix diagonalization.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
static void computeFAndDCosts(const std::vector<Eigen::MatrixXd> &CVec,double &DCost, double &FCost) {
    const size_t N=CVec[0].rows();
    const size_t K=CVec.size();
    double ND=static_cast<double>(N);
    double NProd=(ND*(ND-1));
    
    //Test for convergence.
    FCost=0;
    DCost=0;
    for(size_t k=0;k<K;k++) {
        const Eigen::MatrixXd CCur=CVec[k];
        //FCost holds the cumulative squared Frobenius norm of the off-
        //diagonal elements an d DCost holds the norm of the diagonal
        //elements.
        size_t curIdx=0;
        for(size_t i=0;i<N;i++) {
            for(size_t j=0;j<N;j++) {
                if(i!=j) {
                    FCost+=CCur(curIdx)*CCur(curIdx);
                } else {
                    DCost+=CCur(curIdx)*CCur(curIdx);
                }
                curIdx++;
            }
        }
    }
    
    DCost/=ND;
    FCost/=NProd;
}

/**JOINTMATDIAGZIEHEALG Given K real NXN symmetric matrices in A, try to
*               find a common diagonalization matrix V such that
*               V*A[k]*V.transpose() is diagonal for all k. V is generally
*               not orthogonal. When K>2, it is possible that no common
*               matrix exists, in which case, the sum of the squared
*               Frobenius norms of the matrices is minimized. This function
*               implements equation 17 of [1].
*
*INPUTS: AVec A standard template vector array containing the K NXN
*             Eigen::MatrixXd matrices.
*     maxIter The maximum number of iterations to perform. 1000 is a
*             reasonable number. The algorithm is typically much faster.
*  RelTol, AbsTol Tolerances on the cost function for convergence. A value
*             of std::numeric_limits<double>::epsilon() is reasonable.
*           V Space for the NXN return value.
*        CVec A length K vector array into which V*A[k]*V.transpose() is
*             placed.
*       FCost The double value into which the final cost is placed.
*
*OUTPUTS: The return value is 0 if the algorithm converged and 1 if it did
*         not . The results are put into V, CVec, and FCost.
*
*REFERENCES:
*[1] A. Ziehe, P. Laskov, G. Nolte, and K.-R. Müller, "A fast algorithm for
*    joint diagonalization with non-orthogonal transformations and its
*    application to blind source separation," Journal of Machine Learning,
*    vol. 5, pp. 777-800, Jul. 2004.
*
*June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
int JointMatDiagZieheAlg(const std::vector<Eigen::MatrixXd> &AVec,const size_t maxIter,const double RelTol,const double AbsTol, Eigen::MatrixXd &V, std::vector<Eigen::MatrixXd> &CVec, double &FCost) {
const size_t N=AVec[0].rows();
const size_t K=AVec.size();

//The bound for the Frobenius norm of W. This is an arbitrary value <1.
const double theta=0.9;

Eigen::MatrixXd z(N,N);
Eigen::MatrixXd y(N,N);

//Initial values
Eigen::MatrixXd W=Eigen::MatrixXd::Zero(N,N);
V.setIdentity(N,N);
CVec.resize(K);
for(size_t k=0;k<K;k++) {
    CVec[k]=AVec[k];
}

for(size_t curIter=0;curIter<maxIter;curIter++) {
    //Implementation of y and z for Equation 17.
    z.setZero(N,N);
    y.setZero(N,N);
    for(size_t i=0;i<N;i++) {
        for(size_t j=0;j<N;j++) {
            for(size_t k=0;k<K;k++) {
                z(i,j)=z(i,j)+CVec[k](i,i)*CVec[k](j,j);
                y(i,j)=y(i,j)+CVec[k](j,j)*(CVec[k](i,j)+CVec[k](j,i))/2;
            }
        }
    }

    //Equation 17
    for(size_t i=0;i<N;i++) {
        for(size_t j=(i+1);j<N;j++) {
            const double denom=z(j,j)*z(i,i)-z(i,j)*z(i,j);
            
            W(i,j)=(z(i,j)*y(j,i)-z(i,i)*y(i,j))/denom;
            W(j,i)=(z(i,j)*y(i,j)-z(j,j)*y(j,i))/denom;
        }
    } 

    //Limit the Frobenius norm of W.
    const double WNorm=W.norm();    
    if(WNorm>theta) {
         W=(theta/WNorm)*W;
    }
    
    //Equation 7
    const Eigen::MatrixXd IWSum=W+Eigen::MatrixXd::Identity(N,N);
    V=IWSum*V;
    
    for(size_t k=0;k<K;k++) {
        CVec[k]=IWSum*CVec[k]*IWSum.transpose();
    }
    
    //Test for convergence.
    double DCost;
    computeFAndDCosts(CVec, DCost, FCost);
    if(FCost<=AbsTol||FCost<=RelTol*DCost) {
        return 0;
    }
}
return 1;
}

/**JOINTMATDIAGFROBVOLLALG Given K real NXN symmetric matrices in A, try to
*               find a common diagonalization matrix W such that
*               W*A[k]*W.transpose() is diagonal for all k. W is generally
*               not orthogonal. When K>2, it is possible that no common
*               matrix exists, in which case, the sum of the weighted (by
*               w) squared Frobenius norms of the matrices is minimized.
*               This function implements the algorithm of [1].
*
*INPUTS: AVec A standard template vector array containing the K NXN
*             Eigen::MatrixXd matrices.
*           w A length K weighting array such that the elements sum to 1.
*     maxIter The maximum number of iterations to perform. 1000 is a
*             reasonable number. The algorithm is typically much faster.
*  RelTol, AbsTol Tolerances on the cost function for convergence. A value
*             of std::numeric_limits<double>::epsilon(); is reasonable.
*           W Space for the NXN return value.
*        CVec A length K vector array into which W*A[k]*W.transpose() is
*             placed.
*       FCost The double value into which the final cost is placed.
*
*OUTPUTS: The return value is 0 if the algorithm converged and 1 if it did
*         not. The results are put into W, CVec, and FCost. If NaNs or Inf
*         values were encountered before the eigenvalue decomposition, the
*         return value is 100 and NaNs are put into FCost and W; CVec will
*         not hold anything in particular.
*
*REFERENCES:
*[1] R. Vollgraf and K. Obermayer, "Quadratic optimization for simultaneous
*    matrix diagonalization," IEEE Transactions on Signal Processing, vol.
*    54, no. 9, pp. 3270-3278, Sep. 2006.
*
*June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
int JointMatDiagFrobVollAlg(const std::vector<Eigen::MatrixXd> &AVec,const Eigen::MatrixXd &w, const size_t maxIter,const double RelTol,const double AbsTol, Eigen::MatrixXd &W, std::vector<Eigen::MatrixXd> &CVec, double &FCost) {
    const size_t N=AVec[0].rows();
    const size_t K=AVec.size();

    //We just set C0=eye(N,N) and thus P=eye(N,N) so P'*C0*P=eye(N,N);
    //There does not seem to be a reason to randomize the W matrix, so we
    //just set it to eye(N,N) and thus diag(W*C0*W')=1 for all i.
    W.setIdentity(N,N);
    Eigen::MatrixXd M=Eigen::MatrixXd::Zero(N,N); 
    for(size_t k=0;k<K;k++) {
        //W is the identity matrix, so it can be omitted here.
        const Eigen::MatrixXd M1=AVec[k];
        const Eigen::MatrixXd M1M1p=M1*M1.transpose();
    
        M=M+w(k)*(M1M1p+M1M1p.transpose());
    }
        
    //Initialize.
    CVec.resize(K);
    for(size_t k=0;k<K;k++) {
        CVec[k].setZero(N,N);
    }
    
    for(size_t curIter=0;curIter<maxIter;curIter++) {
        for(size_t i=0;i<N;i++) {
            Eigen::MatrixXd wi=W.col(i);
            for(size_t k=0;k<K;k++) {
                const Eigen::MatrixXd Ak=AVec[k];
                const Eigen::MatrixXd m1=Ak*wi;
                const Eigen::MatrixXd m2=Ak.transpose()*wi;
                M=M-w(k)*(m1*m1.transpose()+m2*m2.transpose());
            }
            
            //Test for non-finite values before calling EigenSolver.
            //Otherwise, EigenSolver won't set the outputs and one can read
            //past the end of the matrices.
            if(anyNonFiniteInMatrix(M)) {
                W=MatrixOfVal<Eigen::MatrixXd,double>(NAN,N,N);
                FCost=NAN;
                return 100;
            }

            //Find the index of the minimum magnitude eigenvalue.
            const Eigen::EigenSolver<Eigen::MatrixXd> eigenM(M);
            Eigen::MatrixXcd d=eigenM.eigenvalues();
            
            size_t minIdx,otherIdx;
            d.cwiseAbs().minCoeff(&minIdx,&otherIdx);
            const Eigen::MatrixXcd V=eigenM.eigenvectors();
            W.col(i)=V.col(minIdx).real();
            wi=W.col(i);

            for(size_t k=0;k<K;k++) {
                const Eigen::MatrixXd Ak=AVec[k];
                const Eigen::MatrixXd m1=Ak*wi;
                const Eigen::MatrixXd m2=Ak.transpose()*wi;
                M=M+w(k)*(m1*m1.transpose()+m2*m2.transpose());
            }
        }
    
        //Compute the LW cost in Equation 4. We only do this every few
        //iterations so that it is faster.
        if(((curIter+1)%25)==0) {//
            double DCost;
            
            for(size_t k=0;k<K;k++) {
                CVec[k]=W.transpose()*AVec[k]*W;
            }
            
            computeFAndDCosts(CVec, DCost, FCost);
            
            if(FCost<=AbsTol||FCost<=RelTol*DCost) {
                W.transposeInPlace();
                return 0;//It converged.
            }
        }
    }

    for(size_t k=0;k<K;k++) {
        CVec[k]=W.transpose()*AVec[k]*W;
    }

    double DCost;
    computeFAndDCosts(CVec, DCost, FCost);
    W.transposeInPlace();
    return 1;  
}

/**APPROXINVMAT1NORM This function approximates the 1-norm of the inverse
 *     of the matrix X without actually performing matrix inversion and can
 *     still give a usuable result with singular and near singular
 *     matrices. v ends up being an approximate null vector. This function
 *     is used in the invCondNumberApprox function and implements the
 *     algorithm of [1].
 *
 *Note that X should not contain any Inf or NaN values, because some
 *versions of Eigen::FullPivLU, which is called in this function, do not
 *allocate space for their return values if that is the case.
 *
 *REFERENCES:
 *[1] J. Higham, Nicholas, "FORTRAN codes for estimating the one-norm
 *   of a real or complex matrix with applications to condition
 *   estimation," ACM Transactions on Mathematical Software, vol. 14, no.
 *   4, pp. 381–396, Dec. 1988.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
static double approxInvMat1Norm(const Eigen::MatrixXd &X, Eigen::MatrixXd &v) {
    const size_t N=X.rows();
    const double ND=static_cast<double>(N);
    
    const Eigen::FullPivLU<Eigen::MatrixXd> LUX(X);
    Eigen::MatrixXd L=Eigen::MatrixXd::Identity(N,N);
    L.triangularView<Eigen::StrictlyLower>()=LUX.matrixLU();
    const Eigen::MatrixXd U=LUX.matrixLU().triangularView<Eigen::Upper>();
    //L must be lower triangular and U upper triangular. The permutations
    //can be ignored.
                
    //Initialize v to a vector full of 1/N.
    v=Eigen::MatrixXd(N,1);
    {
        const double NDInv=1/ND;
        for(size_t k=0;k<N;k++) {
            v(k)=NDInv;
        }
    }

    forwardSubstitutionByCol<Eigen::MatrixXd>(L,v);
    //The result is placed in v.
    backSubstitutionByCol<Eigen::MatrixXd>(U,v); 

    //Get the l1 norm.
    double gammaVal=v.lpNorm<1>();

    if(N==1) {
        v=v.cwiseAbs().eval();
        return gammaVal;
    }

    Eigen::MatrixXd zeta(N,1);
    for(size_t k=0;k<N;k++) {
        //zeta is -1 or 1.
        zeta(k)=2.0*(v(k)>=0)-1.0;
    }
    
    Eigen::MatrixXd x=zeta;
    forwardSubstitutionByCol<Eigen::MatrixXd>(U.transpose(), x);
    //The result is put in x.
    backSubstitutionByCol<Eigen::MatrixXd>(L.transpose(),x); 

    //j is the index of the maximum absolute value of x.
    double maxVal=fabs(x(0));
    size_t j=0;
    for(size_t k=1;k<N;k++) {
        const double curVal=fabs(x(k));
        if(curVal>maxVal) {
            maxVal=curVal;
            j=k;
        }
    }

    size_t k=1;
    const Eigen::MatrixXd I=Eigen::MatrixXd::Identity(N,N);
    Eigen::MatrixXd zetaNew(N,1);
    while(1) {
        const double gammaBar=gammaVal;
        
        v=I.col(j);
        forwardSubstitutionByCol<Eigen::MatrixXd>(L, v);
        //The result is put in v.
        backSubstitutionByCol<Eigen::MatrixXd>(U,v);
        gammaVal=v.lpNorm<1>();
        
        for(size_t i=0;i<N;i++) {
            //zetaNew is -1 or 1.
            zetaNew(i)=2.0*(v(i)>=0)-1.0;
        }
        
        bool allEqual=true;
        for(size_t i=0;i<N;i++) {
            if(zetaNew(i)!=zeta(i)) {
                allEqual=false;
                break;
            }
        }

        if(allEqual||gammaVal<gammaBar) {
            break;
        }
        
        zeta=zetaNew;
        
        x=zeta;
        forwardSubstitutionByCol<Eigen::MatrixXd>(U.transpose(), x);
        //The result is put in x.
        backSubstitutionByCol<Eigen::MatrixXd>(L.transpose(),x); 
        k++;
        
        if(k>4) {
            break;
        }
        
        //This is [~,jNew]=max(abs(x));
        maxVal=fabs(x(0));
        size_t jNew=0;
        for(size_t i=1;i<N;i++) {
            const double curVal=fabs(x(i));
            if(curVal>maxVal) {
                maxVal=curVal;
                jNew=i;
            }
        }

        if(jNew==j) {
            break;
        }
        j=jNew;
    }
    const double ND1=ND-1;

    double signVal=1.0;
    double idx=0;
    for(size_t i=0;i<N;i++) {
        x(i)=signVal*(1.0+idx/ND1);
        signVal=-signVal;
        idx++;
    }

    forwardSubstitutionByCol<Eigen::MatrixXd>(L, x);
    //The result is put in x.
    backSubstitutionByCol<Eigen::MatrixXd>(U,x); 
            
    double testVal=2.0*x.lpNorm<1>()/(3.0*ND);
    if(testVal>gammaVal) {
        v=x;
        gammaVal=testVal;
    }

    if(std::isnan(gammaVal)) {
        gammaVal=std::numeric_limits<double>::infinity();
    }
    return gammaVal;
}

/**INVCONDNUMBERAPPROX Find the approximate inverse condition number of a
 *     real, square matrix using the l1 norm. The output is very similar
 *     to the rcond function built into Matlab, when given a real matrix.
 *     It approximates the inverse condition number without directly
 *     inverting the matrix X by using the approximation of [1].
 *
 *INPUTS: X An NXN real matrix
 *
 *OUTPUTS: The return value is an approximation of the inverse condition
 *         number of X.
 *
 *If X contains any non-finite values, then a NaN is returned.
 *
 *REFERENCES:
 *[1] J. Higham, Nicholas, "FORTRAN codes for estimating the one-norm
 *   of a real or complex matrix with applications to condition
 *   estimation," ACM Transactions on Mathematical Software, vol. 14, no.
 *   4, pp. 381–396, Dec. 1988.
 *
 *June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
double invCondNumberApprox(const Eigen::MatrixXd &X) {
    const size_t N=X.rows();
 
    //Look for NaNs and Inf values in X. If any are found, then return NaN.
    if(anyNonFiniteInMatrix(X)) {
        return NAN;
    }
    
    //Get the 1-matrix norm of X. This is the maximum absolute column sum
    //of the matrix.
    double normX=X.col(0).cwiseAbs().sum();
    for(size_t i=1;i<N;i++) {
        const double curAbsColSum=X.col(i).cwiseAbs().sum();
        normX=std::max(curAbsColSum,normX);
    }

    if(normX==0) {
        return 0;
    }
    
    Eigen::MatrixXd v;
    return 1.0/(normX*approxInvMat1Norm(X, v));
}

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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
