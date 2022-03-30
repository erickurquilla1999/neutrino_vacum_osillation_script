#include <iostream>
using std::cerr;
using std::endl;
#include <complex>
#include <cmath>
#include <fstream>
using std::ofstream;
#include <cstdlib>
#include <string.h>
using namespace std;
#include "density_matrix.h"


int main(){

    //mixing angles radians
    double theta_12=1e-6*M_PI/180;
    double theta_13=48.3*M_PI/180;
    double theta_23=8.61*M_PI/180;
            
    //neutrino mases kg
    double m1=0.049*1.782662e-36;
    double m2=0*1.782662e-36;
    double m3=0*1.782662e-36;
            
    //cp violation angle radians
    double delta_cp=222*M_PI/180;

    //energia J
    double E=10.0e6*1.60218e-19; 

    //size step
    double dt=1.0e-6;

    Density_matrix density;
    std::cout<<real(density.density_matrix[0][0])<<std::endl;    
    for (int i=0;i<100;i++){
        density.evolve_density_matrix(dt, theta_12, theta_13, theta_23, m1, m2, m3, delta_cp, E);
        std::cout<<real(density.density_matrix[0][0])<<std::endl;
    };

    return 0;
};