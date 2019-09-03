#include <iostream>
#include <math.h>
#include <fstream>

void statespace(double* f, double* y);

//A routine for explicit RK4
int main(int argc, char const *argv[]) {
  // std::ofstream datalog;
  double y[6] = {0,M_PI/2,0,9,0,0}, f[6]={0,0,0,0,0,0}, yt[6], yout[6];
  double h=1.0/100.0, hh=h*0.5, h6=h/6.0;
  double k1[6], k2[6], k3[6], k4[6];
  int n=6;
  // datalog.open ("datalogger.csv");
  for (int j = 0; j < 1000; j++){
    // datalog<<y[0]<<","<<y[1]<<","<<y[2]<<","<<y[3]<<","<<y[4]<<","<<y[5]<< "\n";
    std::cout<<y[0]<<","<<y[1]<<","<<y[2]<<","<<y[3]<<","<<y[4]<<","<<y[5]<< "\n";
  	// std::cout<<"yin: "<<y[0]<<" "<<y[1]<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<y[5]<<std::endl;
    statespace(k1, y);
    // std::cout<<" "<<std::endl;
    for (int i = 0; i < n; i++) yt[i]=y[i]+hh*k1[i];
    // std::cout<<hh<<std::endl;
    // std::cout<<"k1: "<<k1[0]<<" "<<k1[1]<<" "<<k1[2]<<" "<<k1[3]<<" "<<k1[4]<<" "<<k1[5]<<std::endl;
    // std::cout<<"yt: "<<yt[0]<<" "<<yt[1]<<" "<<yt[2]<<" "<<yt[3]<<" "<<yt[4]<<" "<<yt[5]<<std::endl;
    statespace(k2, yt);
    // std::cout<<" "<<std::endl;
    for (int i = 0; i < n; i++) yt[i]=y[i]+hh*k2[i];
    // std::cout<<"k2: "<<k2[0]<<" "<<k2[1]<<" "<<k2[2]<<" "<<k2[3]<<" "<<k2[4]<<" "<<k2[5]<<std::endl;
    // std::cout<<"yt: "<<yt[0]<<" "<<yt[1]<<" "<<yt[2]<<" "<<yt[3]<<" "<<yt[4]<<" "<<yt[5]<<std::endl;
    statespace(k3, yt);
    // std::cout<<" "<<std::endl;
    for (int i = 0; i < n; i++) {
      yt[i]=y[i]+h*k3[i];
      k3[i]+=k2[i];
    };
    // std::cout<<"k3: "<<k3[0]<<" "<<k3[1]<<" "<<k3[2]<<" "<<k3[3]<<" "<<k3[4]<<" "<<k3[5]<<std::endl;
    // std::cout<<"yt: "<<yt[0]<<" "<<yt[1]<<" "<<yt[2]<<" "<<yt[3]<<" "<<yt[4]<<" "<<yt[5]<<std::endl;
    statespace(k4, yt);
    // std::cout<<" "<<std::endl;
    for (int i = 0; i < n; i++) {
      yout[i]=y[i]+h6*(k1[i]+2*k3[i]+k4[i]);
  	  y[i]=yout[i];
    }
    // std::cout<<"k4: "<<k4[0]<<" "<<k4[1]<<" "<<k4[2]<<" "<<k4[3]<<" "<<k4[4]<<" "<<k4[5]<<std::endl;
    // std::cout<<" "<<std::endl;
  	// std::cout<<"yout: "<<yout[0]<<" "<<yout[1]<<" "<<yout[2]<<" "<<yout[3]<<" "<<yout[4]<<" "<<yout[5]<<std::endl;
    // std::cout<<yout[1], yout[2], yout[3], yout[4]<<std::endl;
  }
  // datalog.close();
  return 0;
}

//y[] contains {theta, phi, psi, dtheta, dphi, dpsi} respectively
//state-space returns
void statespace(double* f, double* y){
double cp = cos(y[1]), sp = sin(y[1]), cs = cos(y[2]), ss = sin(y[2]), dt = y[3],\
        dp = y[4], ds = y[5];
  f[0] = dt;
  f[1] = dp;
  f[2] = ds;
  f[3] = (-5*cs*dt*(ds + cp*dt)*ss)/3. - (dp*pow(sp,-1)*(ds*(1 + 5*pow(cs,2) - 5*pow(ss,2)) + cp*dt*(13 + 5*pow(cs,2) - 5*pow(ss,2))))/6.;
  f[4] = (ds*(10*cs*dp*ss + dt*sp*(1 - 5*pow(cs,2) + 5*pow(ss,2))) + cp*dt*(10*cs*dp*ss + dt*sp*(7 - 5*pow(cs,2) + 5*pow(ss,2))))/6.;
  f[5] = (pow(sp,-1)*(-60*cs*sp*ss*pow(dp,2) + 10*cs*sp*ss*pow(dt,2)*(5 + pow(sp,2)) + 2*cp*ds*(dp + 10*cs*dt*sp*ss + 5*dp*pow(cs,2) - 5*dp*pow(ss,2)) + dp*dt*(19 + 5*pow(cs,2)*(7 + 5*pow(sp,2)) - 35*pow(ss,2) - pow(sp,2)*(7 + 25*pow(ss,2))) + dt*pow(cp,2)*(-30*cs*dt*sp*ss + dp*(7 - 25*pow(cs,2) + 25*pow(ss,2)))))/12.;
}
