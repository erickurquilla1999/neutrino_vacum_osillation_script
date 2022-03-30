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



Density_matrix::Density_matrix()
{
    num_steps=0;
    density_matrix[0][0]={1,0};

    //save initial conditions density matrix
    ofstream outdata;
    outdata.open("/home/centroescolarjuanabarrera/tesis/v_o_s/output/step_0.dat"); // opens the file
        if( !outdata ) 
        {
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }
    string outvalue;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            double re=real(density_matrix[i][j]);
            double im=imag(density_matrix[i][j]);
            outvalue=to_string(re)+" "+to_string(im)+" 0.0";
            outdata << outvalue << endl;
        };  
    };  
    outdata.close(); 
};





double Density_matrix::evolve_density_matrix(double dt,double theta_12,double theta_13,double theta_23,double m1,double m2,double m3,double delta_cp,double E)
{

    //light velocity
    //double c=1.0;
    double c=299792458;

    //hbar 
    //double hbar=1/(2*M_PI);
    double hbar=1.054571817e-34;

    // sin and cos for U matix
    double s12=sin(theta_12);
    double c12=cos(theta_12);
    double s13=sin(theta_13);
    double c13=cos(theta_13);
    double s23=sin(theta_23);
    double c23=cos(theta_23);
    double sdcp=sin(delta_cp);
    double cdcp=cos(delta_cp);

    //obtain mixing matrix U
    std::complex<double> mixing_matrix[3][3]={
    { {c12*c13,0} , {s12*c13,0} , {s13*cdcp,-s13*sdcp} } ,
    { {-s12*c23-c12*s13*s23*cdcp,-c12*s13*s23*sdcp} , {c12*c23-s12*s13*s23*cdcp,-s12*s13*s23*sdcp} , {c13*s23,0} },
    { {s12*s23-c12*s13*c23*cdcp,-c12*s13*c23*sdcp} , {-c12*s23-s12*s13*c23*cdcp,-s12*s13*c23*sdcp} , {c13*c23,0} }
    };

    //obtain mixin matrix U daguer
    std::complex<double> mixing_matrix_daguer[3][3];
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            mixing_matrix_daguer[j][i]=std::conj(mixing_matrix[i][j]);
        };
    };

    //mass matrix
    std::complex<double> mass_matrix[3][3];
    mass_matrix[0][0]={m1*m1*c*c*c*c/(2*E),0};
    mass_matrix[1][1]={m2*m2*c*c*c*c/(2*E),0};
    mass_matrix[2][2]={m3*m3*c*c*c*c/(2*E),0};

    // obtain hamiltonian
    std::complex<double> hamiltonian[3][3];
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            hamiltonian[i][j]=0;
            for (int a=0;a<3;a++){
                for (int b=0;b<3;b++){
                    hamiltonian[i][j]=hamiltonian[i][j]+mixing_matrix[i][a]*mass_matrix[a][b]*mixing_matrix_daguer[b][j];
                    };  
            };
        };
    };

    //######################################
    //doing a step with RK4
    //######################################

    std::complex<double> k1[3][3];
    std::complex<double> k2[3][3];
    std::complex<double> k3[3][3];
    std::complex<double> k4[3][3];

    complex<double> minus_i_over_hbar={0,-1/hbar};
            
    //found k1
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            k1[i][j]=0;
            for (int a=0;a<3;a++){
                k1[i][j]=k1[i][j]+hamiltonian[i][a]*density_matrix[a][j]-density_matrix[i][a]*hamiltonian[a][j];
            };
            k1[i][j]=minus_i_over_hbar*k1[i][j];
        };
    };    

    //found k2
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            k2[i][j]=0;
            for (int a=0;a<3;a++){
                k2[i][j]=k2[i][j]+hamiltonian[i][a]*(density_matrix[a][j]+0.5*k1[a][j]*dt)-(density_matrix[i][a]+0.5*k1[i][a]*dt)*hamiltonian[a][j];
            };
            k2[i][j]=minus_i_over_hbar*k2[i][j];
        };     
    };

    //found k3
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            k3[i][j]=0;
            for (int a=0;a<3;a++){
                k3[i][j]=k3[i][j]+hamiltonian[i][a]*(density_matrix[a][j]+0.5*k2[a][j]*dt)-(density_matrix[i][a]+0.5*k2[i][a]*dt)*hamiltonian[a][j];
            };
            k3[i][j]=minus_i_over_hbar*k3[i][j];
        };
    };  

    //found k4
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            k4[i][j]=0;
            for (int a=0;a<3;a++){
                k4[i][j]=k4[i][j]+hamiltonian[i][a]*(density_matrix[a][j]+k3[a][j]*dt)-(density_matrix[i][a]+k3[i][a]*dt)*hamiltonian[a][j];
            };
            k4[i][j]=minus_i_over_hbar*k4[i][j];
        };
    };  

    //found next step density matrix 
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            density_matrix[i][j]=density_matrix[i][j]+(1.0/6.0)*dt*(k1[i][j]+2.0*k2[i][j]+2.0*k3[i][j]+k4[i][j]);
        };
    };  

    //sum one to the count step variable
    num_steps=num_steps+1;

    //save data
    ofstream outdata;
    outdata.open("/home/centroescolarjuanabarrera/tesis/v_o_s/output/step_"+to_string(num_steps)+".dat"); // opens the file
        if( !outdata ) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
                exit(1);
        }

    string outvalue;
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            double re=real(density_matrix[i][j]);
            double im=imag(density_matrix[i][j]);
            double time=num_steps*dt;
            outvalue=to_string(re)+" "+to_string(im)+" "+to_string(time);
            outdata << outvalue << endl;
        };  
    };  
    outdata.close();

    return 0,0;

};
