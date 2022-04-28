/*
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;




int main()
{
    
  complex<double> z0=1+1i;
  cout<<"z0= "<< z0 << endl;


 
  complex<double> z1 = 1i * 1i;     // imaginary unit squared
  cout << "i * i = " << z1 << '\n';
 
  complex<double> z2 = pow(1i, 2); // imaginary unit squared
  cout << "pow(i, 2) = " << z2 << '\n';
 


  const double PI = M_PI; // or std::numbers::pi in C++20
  complex<double> z3 = exp(1i * PI); // Euler's formula
  cout << "exp(i * pi) = " << z3 << '\n';
 
  complex<double> z4 = 1. + 2i, z5 = 1. - 2i; // conjugates
  cout << "(1+2i)*(1-2i) = " << z4*z5 << '\n';

    return 0;
}


// complex_op_as.cpp
// compile with: /EHsc
#include <complex>
#include <iostream>

int main( )
{
   using namespace std;
   double pi = 3.14159265359;

   // Example of the first member function
   // type complex<double> assigned to type complex<double>
   complex<double> cl1( 3.0 , 4.0 );
   complex<double> cr1( 2.0 , -1.0 );
   cout << "The left-side complex number is cl1 = " << cl1 << endl;
   cout << "The right-side complex number is cr1 = " << cr1 << endl;

   cl1  = cr1;
   cout << "The complex number cr1 assigned to the complex number cl1 is:"
        << "\ncl1 = cr1 = " << cl1 << endl;

   // Example of the second member function
   // type double assigned to type complex<double>
   complex<double> cl2( -2 , 4 );
   double cr2 =5.0;
   cout << "The left-side complex number is cl2 = " << cl2 << endl;
   cout << "The right-side complex number is cr2 = " << cr2 << endl;

   cl2 = cr2;
   cout << "The complex number cr2 assigned to the complex number cl2 is:"
        << "\ncl2 = cr2 = " << cl2 << endl;

   cl2 = complex<double>(3.0, 4.0);
   cout << "The complex number (3, 4) assigned to the complex number cl2 is:"
        << "\ncl2 = " << cl2 << endl;
}
*/

#include <iostream>
#include <complex>
using namespace std;
int main() {
   complex<double> complexnumber(4, 3);
   cout<<"Square of (4+ 3i) = "<<pow(complexnumber, 2)<<endl;
   
complex<double> cl2( -2 , 4 );
   double cr2 =5.0;
   cout << "The left-side complex number is cl2 = " << cl2 << endl;
   cout << "The right-side complex number is cr2 = " << cr2 << endl;

   cl2 = cr2;
   cout << "The complex number cr2 assigned to the complex number cl2 is:"
        << "\ncl2 = cr2 = " << cl2 << endl;

   cl2 = complex<double>(3.0, 4.0);
   cout << "The complex number (3, 4) assigned to the complex number cl2 is:"
        << "\ncl2 = " << cl2 << endl;


   return 0;
}