///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  double k_linear = p->b*p->d;
  double k_quad = p->c*p->d*p->d;
  return -k_linear*y[1]/p->m - k_quad * sqrt(y[1]*y[1] + y[3]*y[3]+y[5]*y[5]) * y[1] / p->m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

double f_rj(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[3];
}

double f_vj(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  double k_linear = p->b*p->d;
  double k_quad = p->c*p->d*p->d;
  return -k_linear*y[3]/p->m - k_quad * sqrt(y[1]*y[1] + y[3]*y[3]+y[5]*y[5]) * y[3] / p->m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

double f_rk(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[5];
}

double f_vk(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  double k_linear = p->b*p->d;
  double k_quad = p->c*p->d*p->d;
  return -k_linear*y[5]/p->m - k_quad * sqrt(y[1]*y[1] + y[3]*y[3]+y[5]*y[5]) * y[5] / p->m - p->g;
  // return -g;    // if no air constant acceleration along -j direction: F/m = -g
}
 
double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  if (y[0]>18.5) return 1; 
  return 0;  // continue calculation
}

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.0075;   
  pars.b=1.6e-4;  
  pars.c=0.25;
  void *p_par = (void*) &pars;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  double vPitch = 45.2001;
  while ((c = getopt (argc, argv, "x:z:t:v:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'v':
      vPitch = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
      break;
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


     // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here

  vector<pfunc_t> v_fun(6);  
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  v_fun[4]=f_rk;
  v_fun[5]=f_vk;

  vector<double> y(6);
  // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
  y[0]=0;   // init position on i-axis
  y[1]=vPitch*cos(theta0*3.14159/180);  // init velocity along i axis
  y[2]=0;   // repeat for j-axis
  y[3]=0;
  y[4]=z0;
  y[5]=vPitch*sin(theta0*3.14159/180);
  //cout << "Vinit: " << vPitch << " m/s" << endl;
  //cout << "Angle: " << theta0 << " deg" << endl;
  //cout << "(vx,vx) " << y[1] << " , "  <<  y[5] << " m/s" << endl;

  double x=0;           // t0
  double xmax=5;  // tmax
  int nsteps=100000000;
  // fixed step size algorithm
  auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

  //int N = tgN[0].GetN();
  //double t_last,ri_last,rk_last;
  //int i=1;
  //double curxval=19.0;
  /*
  for (int i=1;i<6;i++){
    tgN[0].GetPoint(N-i,t_last,ri_last);
    tgN[4].GetPoint(N-i,t_last,rk_last);
    cout<<"x pos: "<<ri_last<<endl;
    cout<<"z pos: "<<rk_last<<endl;
    cout<<"----------"<<endl;
    //curxval=ri_last;
    
    }*/
  
  
  //cout<<"x pos: "<<y[0]<<endl;
  //cout<<"z pos: "<<y[4]<<endl;

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

