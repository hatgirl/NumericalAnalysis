// Runge Kutta 4 solver specifically for the 2nd order ODE given
// in the 3rd Numerical Analysis project, relating pressure gradients
// of a star to radius from the center of the star.
//
// @author Jessica Sorrell
// @date 29-Nov-2014

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

// order of x
float N = 1;

// returns -x^n
double xNFunc( double t, double dxdt, 
		    double dvdt){
  return ((2/t)*dxdt + dvdt);
}



// returns x
double xFunc( double t, double dxdt, 
		  double dvdt){
  if( t > pow(10,-7)){
    return (-pow(xNFunc( t, dxdt, dvdt), 1/N));
  }
  else {
    return (1 - pow(t, 2)/6);
  }
}


// returns double derivative of x with respect to squiggle
double dvdtFunct (double t, double x, double dxdt){
  return (-(pow(x, N)) - (2 / t)*dxdt);
}



//repopulates the nextSteps array with next values for x and dx
void steps( double timeStep, double t, double x,
		 double dxdt, double dvdt, 
		 double* nextSteps ){
  
  double dx1 = timeStep*dxdt;
  double dxdt1 = 
    timeStep*dvdtFunct(t, x, dxdt );
  
  double dx2 = timeStep*(dxdt + dxdt1 / 2);
  double dxdt2 = 
    timeStep*dvdtFunct( (t + timeStep/2), 
			      (x + dx1/2 ),
			      (dxdt + dxdt1/2));
  
  double dx3 = timeStep*(dxdt + dxdt2 / 2);
  double dxdt3 = 
    timeStep*dvdtFunct( (t + timeStep/2), 
			      (x + dx2/2 ),
			      (dxdt + dxdt1/2));
  

  double dx4 = timeStep*(dxdt + dxdt3);
  double dxdt4 = 
    timeStep*dvdtFunct( (t + timeStep), 
			      (x + dx3 ),
			      (dxdt + dxdt1));

  double dx = (dx1 + 2*dx2 + 2*dx3 + dx4)/6;
  double dv = 
    (dxdt1 + 2*dxdt2 + 2*dxdt3 + dxdt4)/6;

  nextSteps[0] = dx;
  nextSteps[1] = dv;

}


int main( int argc, char* argv[]){

  // open file to dump output
  std::ofstream out;
  std::streambuf  *coutbuf;
  out.open("rk4_test.txt");

  coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); 

  // user passes step size
  double timeStep = atof(argv[1]);

  // pressure, pushing down on me, pushing down on you
  double pressure = 0;

 
  double t = 0.00001;
  double dxdt = -t/3;
  double x = 1 - pow(t,2)/6;
  double dvdt = dvdtFunct(t, x, dxdt);

  double step [2] = {0.0, 0.0};

  // iterate until you've reached the surface of the star  
  while( x > 0.0000001 ){

    steps( timeStep, t, x, dxdt, dvdt, step);
    x = x + step[0];
    dxdt = dxdt + step[1];
    t = t + timeStep;
    dvdt = dvdtFunct( t, x, dxdt );

    pressure += x*timeStep;

    // std::cout << "X = \t " << std::setprecision(15) 
    // 	      << x << std::endl;
    
    // std::cout << "X dot = \t " << std::setprecision(15) 
    // 	      << dxdt << std::endl;

    // print t, x, x^n
    std::cout  << std::setprecision(15) << std::setw(15) << t << "\t " <<
      x << "\t" << pow(x, N) << std::endl;
    
    // std::cout << "X doubledot = \t " << std::setprecision(15)
    // 	      << xDoubleDot << std::endl;

    // std::cout << "Pressure = \t " << std::setprecision(15)
    // 	      << xDoubleDot << std::endl;

    
    
  } 

  //Restore the old buffer
  std::cout.rdbuf(coutbuf);
  
  //Close the output file
  out.close();
 
  return 1;
  
}
