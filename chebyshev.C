// Chebyshev coefficient finder
//
// Jessica Sorrell
// GPL, best believe
// 
// Pass it the number of coefficients you need
//
//
// These headers are so standard one tends not to even think about them
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

double func( double x )
{
return cos(2*M_PI*x);
//return log(1+x);
}


int main( int argc, char* argv[] )
{

//number of Chebyshev modes for the expansion
int num_modes = atoi(argv[1]);


//Chebyshev nodes
double* nodes = new double[num_modes];

for( int i = 0; i < num_modes; i++)
  {
nodes[i] = cos( M_PI*(i + 0.5)/num_modes );
}


//our a_ns for all the Chebyshev modes
double* coeffs = new double[num_modes];


//calculate the a_0 coefficient 
for (int j = 0; j < num_modes; j++ )
  {
coeffs[0] += func( nodes[j] );
}
coeffs[0] /= num_modes;

//print the coefficients to stdout
std::cout << "a_0: \t " << std::setprecision(15) << coeffs[0] << std::endl;


//calculate the a_i coefficients
for (int i = 1; i < num_modes; i++ )
  {

coeffs[i] = 0;


//sum from 0 to n-1 the function at the jth node times the ith 
//Chebyshev interpolating polynomial evaluated at the jth node
for (int j = 0; j < num_modes; j++ )
  {
coeffs[i] += (func( nodes[j] ) * cos(i * acos ( nodes[j] )));
}
coeffs[i] /= 0.5*num_modes;

std::cout << "a_" << i << ": \t " << std::setprecision(15) << coeffs[i] << std::endl;
}


//print two sets of data. one to plot the function. 
//one to plot the fit
  std::ofstream out;
  std::streambuf  *coutbuf;
  out.open("chebyshev.txt");

  coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); 

double fitval;
double val;

for (double x = -1; x < 1; x += 0.02)
  {

fitval = 0.0;
for (int j = 0; j < num_modes; j++ ){

//calculate the values of the function from the chebyshev interp
fitval += coeffs[j]*cos(j*acos(x));
}

//calculate the actual values of the function for comparison
val = func(x);

//print x values, real values, and interp values to a file
std::cout << "x: "<< std::setprecision(12) << std::fixed <<  x <<"  \t val: "<< val << 
  "\t fitval: "<< fitval << std::endl;
}

for (int i = 0; i < num_modes; i ++ )
  {
//print the coefficients to a file
std::cout << "a_" << i << ": " << std::setprecision(12) << coeffs[i] << std::endl;
}

  //Restore the old buffer
  std::cout.rdbuf(coutbuf);
  
  //Close the output file
  out.close();

return 1;
}
