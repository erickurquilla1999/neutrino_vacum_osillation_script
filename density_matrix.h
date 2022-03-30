#pragma once
#include <iostream>
using std::cerr;
using std::endl;
#include <complex>
#include <cmath>
#include <fstream>
using std::ofstream;
#include <cstdlib> 
#include <string.h>

//improve_01_change
//improve_01_change
//improve_01_change
//improve_01_change
//improve_01_change

class Density_matrix{
    public:

        //definition of density matrix as a complex array
        std::complex<double> density_matrix[3][3]; 
        //numbers of steps doing 
        int num_steps;
        //initial conditions function
        Density_matrix();
        //function for evolve the density matrix a step with in dt
        double evolve_density_matrix(double dt,double theta_12,double theta_13,double theta_23,double m1,double m2,double m3,double delta_cp,double E);
    };
