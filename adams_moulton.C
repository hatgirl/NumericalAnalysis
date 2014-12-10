// Adams-Moulton 4-step solver for the 1st order ODE
// given in the 4th Numerical Analysis project
//
// @author Jessica Sorrell
// @date 10-Dec-2014

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

// returns double derivative of x with respect to squiggle
double dxdtFunct (double t, double x){

  return t*exp(3.0*t) - 2.0*x;
}


// use RK4 to populate the history array
double startSteps( double timeStep, double t, double x,
		   double dxdt){
  
  double k1 = dxdtFunct(t, x);
  
  double k2 = dxdtFunct(t + timeStep/2.0, x + k1*timeStep/2.0);
   
  double k3 = dxdtFunct(t + timeStep/2.0, x + k2*timeStep/2.0);
   
  double k4 = dxdtFunct(t + timeStep, x + k3*timeStep);
   
  double dx = (k1 + 2*k2 + 2*k3 + k4)/6.0;
  return dx;
}

// Adams-Bashforth 5-step
void predict( double timeStep, double t, double* history ){

  double term1 = 1901*dxdtFunct( t - timeStep, history[4]);
  double term2 = 2774*dxdtFunct( t - 2*timeStep, history[3]);
  double term3 = 2616*dxdtFunct( t - 3*timeStep, history[2]);
  double term4 = 1274*dxdtFunct( t - 4*timeStep, history[1]);
  double term5 = 251*dxdtFunct( t - 5*timeStep, history[0]);
  
  double x5 = history[4] + 
    (timeStep/720.)*(term1 - term2 + term3 - term4 + term5);
  
  for (int i = 0; i < 4; i++ ){
    history[i] = history[i+1];
  }
  
  history[4] = x5;
  std::cout<< x5 << " \t time: \t" << t << std::endl;

}

// Adams-Moulton 4-step
void steps( double timeStep, double t, double* history ){

  predict( timeStep, t, history );
  double term1 = 251*dxdtFunct( t, history[4] );
  double term2 = 646*dxdtFunct( t - timeStep, history[3]);
  double term3 = 264*dxdtFunct(t - 2*timeStep, history[2]);
  double term4 = 106*dxdtFunct(t - 3*timeStep, history[1]);
  double term5 = 19*dxdtFunct(t - 4*timeStep, history[0]);

  double x5 = history[3] +
    (timeStep/720.)*(term1 + term2 - term3 + term4 - term5);

  history[4] = x5;
  std::cout<< x5 << " \t time: \t" << t << std::endl;
}

int main( int argc, char* argv[]){

  // open file to dump output
  std::ofstream out;
  std::streambuf  *coutbuf;
  out.open("am4.txt");

  coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); 

  // user passes step size
  double timeStep = atof(argv[1]);
 
  double t = 0.25;
  double x = 0.00;
  double dxdt = dxdtFunct(t, x); 

  //  populate history with actual values;
  double history[5] = 
    {0.0 , 0.00133847, 0.00575205, 0.0139496, 0.0268128};
  
  //let RK4 populate history
  
  //double history[5] = {x, 0.0, 0.0, 0.0, 0.0};
  //for( int i = 0; i < 4; i++){
    
  //  x += startSteps( timeStep, t, x, dxdt);
  //  history[i+1] = x;
  //  t += timeStep;
  //  dxdt = dxdtFunct( t, x);
  // }
  
  
  while (2.05 - t > timeStep/2.0){
    steps( timeStep, t, history );
    t += timeStep;
  }


  std::cout << std::setprecision(15) << std::setw(15) << t 
	     << "\t " << history[4] << std::endl;
  
  //Restore the old buffer
  std::cout.rdbuf(coutbuf);
  
  //Close the output file
  out.close();
 
  return 1;
  
}
