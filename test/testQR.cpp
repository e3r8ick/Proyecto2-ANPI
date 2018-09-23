#include <boost/test/unit_test.hpp>

#include "QRdecomp.hpp"
#include "Pivot.hpp"
#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>


/*
descomposicion:
	1. probar Q*QT= I
	2. R triangular derecha
	3. QR =A
*/

namespace anpi {
  namespace test {
  	template<typename T>
  	void qrTest(){
  		T eps=sqrt(std::numeric_limits<T>::epsilon());
  		//test matrix
  		anpi::Matrix<T> A = { {-1,-2,1},
                          { 2, 0,1},
                          {-1,-1,0}
                          };
       	int r,c;
       	r=A.rows();
       	c=A.cols();
			anpi::Matrix<T> Q=anpi::Matrix<T>(r,c,0.0);
  		anpi::Matrix<T> R=anpi::Matrix<T>(r,c,0.0);
  		//apply qr decomposition
  		anpi::qr(A,Q,R);

  		//create Identity Matrix
  		anpi::Matrix<T> I=anpi::Matrix<T>(r,c,0.0);
  		for(int i=0;i<r;i++){
    		for(int j=0;j<c;j++){
      		I(i,j)=0;
      		I(i,i)=1;
    		}
  		}

  		anpi::Matrix<T> Qt=anpi::Matrix<T>(r,c,0.0);
  		//transpose Q
  		for(int i=0;i<r;i++){
    		for(int j=0;j<c;j++){
      			Qt(j,i)=Q(i,j);
    		}
  		}
  		anpi::Matrix<T> tmp=anpi::Matrix<T>(r,c,0.0);
  		tmp=Q*Qt; 		
  		//Test Q*Qt=I

  		for(int i =0;i<r;++i){
  			for(int j=0;j<c;++j){
  				BOOST_CHECK(std::abs(I(i,j)-tmp(i,j))<eps);
  			}

  		}

      //Test R is upper matrix
  		for(int j=0;j<c;j++){
    		for(int i=j+1;i<r;i++){
    			BOOST_CHECK(std::abs(R(i,j))<eps);
      			
    		}
  		}		

      //test A== Qt*R
  		tmp=Qt*R;
  		for(int i =0;i<r;++i){
  			for(int j=0;j<c;++j){
  				BOOST_CHECK(std::abs(A(i,j)-tmp(i,j))<eps);
  			}

  		}

      //test Solve { {-1,-2,1}, { 2, 0,1}, {-1,-1,0} }*{x,y,z}={1,0,-1}
      std::vector<T> b={1,0,-1};
      std::vector<T> x={0.,0.,0.};
      anpi::solveQR(A,x,b);
      std::vector<T> xRef={-3,4,6};
      
      for(int i =0;i<r;++i){
        BOOST_CHECK(std::abs(x[i]-xRef[i])<eps);
      }

      
  	}
  }
}



  BOOST_AUTO_TEST_SUITE( QR )

BOOST_AUTO_TEST_CASE(QRdecomp) 
{
  anpi::test::qrTest<float>();
  anpi::test::qrTest<double>();
}



BOOST_AUTO_TEST_SUITE_END()