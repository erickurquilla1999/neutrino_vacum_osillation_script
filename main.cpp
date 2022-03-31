#include <iostream>
#include <complex>
#include "density_matrix.h"

int main(){

    //mixing angles radians
    double theta_12=33.82*M_PI/180;
    double theta_13=8.61*M_PI/180;
    double theta_23=48.3*M_PI/180;
    
    //neutrino mases eV
    double m1=1.0;
    double m2=10.0;
    double m3=10.1;
            
    //cp violation angle radians
    double delta_cp=282*M_PI/180;

    //energia eV
    double E=10.0e6; 

    //time size step natural units
    double dt=7502.310814542824;

    //number of step to be computed
    int number_of_steps=10000;
    
    //define an object of Density_matrix class
    Density_matrix density_object;

    //run numb_of_steps times the evolve_density_matrix() functions
    for (int i=0;i<number_of_steps;i++)
    {
        density_object.evolve_density_matrix(dt, theta_12, theta_13, theta_23, m1, m2, m3, delta_cp, E);
    };

    return 0;

};