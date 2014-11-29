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

// order of theta
float N = 2.0;

// returns -theta^n
double thetaNFunct( double radius, double thetaDot, 
		    double doubleThetaDot){
  return ((2/radius)*thetaDot + doubleThetaDot);
}



// returns theta
double thetaFunc( double radius, double thetaDot, 
		  double doubleThetaDot){
  if( radius > pow(10,-7)){
    return (-pow(thetaNFunct( radius, thetaDot, doubleThetaDot), 1/N));
  }
  else {
    return (1 - pow(radius, 2)/6);
  }
}


// returns double derivative of theta with respect to squiggle
double dThetaDotdRadius (double radius, double theta, double thetaDot){
  return (-(pow(theta, N)) - (2 / radius)*thetaDot);
}



//repopulates the nextSteps array with next values for theta and dtheta
void steps( double radStep, double radius, double theta,
		 double thetaDot, double doubleThetaDot, 
		 double* nextSteps ){
  
  double dTheta1 = radStep*thetaDot;
  double dThetaDot1 = 
    radStep*dThetaDotdRadius(radius, theta, thetaDot );
  
  double dTheta2 = radStep*(thetaDot + dThetaDot1 / 2);
  double dThetaDot2 = 
    radStep*dThetaDotdRadius( (radius + radStep/2), 
			      (theta + dTheta1/2 ),
			      (thetaDot + dThetaDot1/2));
  
  double dTheta3 = radStep*(thetaDot + dThetaDot2 / 2);
  double dThetaDot3 = 
    radStep*dThetaDotdRadius( (radius + radStep/2), 
			      (theta + dTheta2/2 ),
			      (thetaDot + dThetaDot1/2));
  

  double dTheta4 = radStep*(thetaDot + dThetaDot3);
  double dThetaDot4 = 
    radStep*dThetaDotdRadius( (radius + radStep), 
			      (theta + dTheta3 ),
			      (thetaDot + dThetaDot1));

  double dT = (dTheta1 + 2*dTheta2 + 2*dTheta3 + dTheta4)/6;
  double dTdot = 
    (dThetaDot1 + 2*dThetaDot2 + 2*dThetaDot3 + dThetaDot4)/6;

  nextSteps[0] = dT;
  nextSteps[1] = dTdot;

}


int main( int argc, char* argv[]){

  // open file to dump output
  std::ofstream out;
  std::streambuf  *coutbuf;
  out.open("rk4.txt");

  coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); 

  // user passes step size
  double radStep = atof(argv[1]);

  // pressure, pushing down on me, pushing down on you
  double pressure = 0;
  double theta = 1;
  double thetaDot = 0;
  double radius = 0.00001;
  double thetaDoubleDot = dThetaDotdRadius(radius, theta, thetaDot);

  double step [2] = {0.0, 0.0};

  // iterate until you've reached the surface of the star  
  while( theta > 0.0000001 ){

    steps( radStep, radius, theta, thetaDot, thetaDoubleDot, step);
    theta = theta + step[0];
    thetaDot = thetaDot + step[1];
    radius = radius + radStep;
    thetaDoubleDot = dThetaDotdRadius( radius, theta, thetaDot );

    pressure += theta*radStep;

    // std::cout << "Theta = \t " << std::setprecision(15) 
    // 	      << theta << std::endl;
    
    // std::cout << "Theta dot = \t " << std::setprecision(15) 
    // 	      << thetaDot << std::endl;

    // print radius, theta, theta^n
    std::cout  << std::setprecision(15) << std::setw(15) << radius << "\t " <<
      theta << "\t" << pow(theta, N) << std::endl;
    
    // std::cout << "Theta doubledot = \t " << std::setprecision(15)
    // 	      << thetaDoubleDot << std::endl;

    // std::cout << "Pressure = \t " << std::setprecision(15)
    // 	      << thetaDoubleDot << std::endl;

    
    
  } 

  //Restore the old buffer
  std::cout.rdbuf(coutbuf);
  
  //Close the output file
  out.close();
 
  return 1;
  
}
