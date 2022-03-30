#include <iostream>
#include <complex>
#include "density_matrix.h"
<<<<<<< HEAD
//master
//master
//master
=======
//improve_01
//improve_01
//improve_01
>>>>>>> improve_01
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

    //time size step
    double dt=1.0e-6;

    //number of step to be compute
    int number_of_steps=100;
    
    //define an object of Density_matrix class
    Density_matrix density_object;

    //run numb_of_steps times the evolve_density_matrix() functions
    for (int i=0;i<number_of_steps;i++)
    {
        density_object.evolve_density_matrix(dt, theta_12, theta_13, theta_23, m1, m2, m3, delta_cp, E);
    };

    return 0;

    //print final message
    std::cout<<"code finish: result can be found in de output directory"<<std::endl;


};