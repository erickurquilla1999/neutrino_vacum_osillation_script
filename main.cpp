#include <iostream>
using std::cerr;
using std::endl;
#include <complex>
#include <cmath>
#include <fstream>
using std::ofstream;
#include <cstdlib>
#include <string>
#include <string.h>
using namespace std;


class Density_matrix{

    public:

        //definition of density matrix complex array
        std::complex<double> density_matrix[3][3];
        int num_steps=0;

        //initial conditions
        Density_matrix(){
            density_matrix[0][0]={1,0};

            //save data
            ofstream outdata;
            outdata.open("/home/centroescolarjuanabarrera/tesis/v_o_s/output/step_0.dat"); // opens the file
                if( !outdata ) { // file couldn't be opened
                    cerr << "Error: file could not be opened" << endl;
                      exit(1);
                }
            string outvalue;
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    double re=real(density_matrix[i][j]);
                    double im=imag(density_matrix[i][j]);
                    outvalue=to_string(re)+" "+to_string(im)+" 0.0";
                    outdata << outvalue << endl;
                };  
            };  
            outdata.close();
        };

        void print_matrix(std::complex<double> matrix[3][3]){
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    std::cout<<matrix[i][j]<<" ";
                };
                std::cout<<""<<std::endl;
            };
        };
        double evolve_density_matrix(double dt){

            //mixing angles radians
            double theta_12=1e-6*M_PI/180;
            double theta_13=48.3*M_PI/180;
            double theta_23=8.61*M_PI/180;

            //cp violation angle radians
            double delta_cp=222*M_PI/180;

            //energia J
            double E=10.0e6*1.60218e-19; 

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

            //neutrino mases kg
            double m1=0.049*1.782662e-36;
            double m2=0*1.782662e-36;
            double m3=0*1.782662e-36;
            
            //obtain mixing matrix
            std::complex<double> mixing_matrix[3][3]={
            { {c12*c13,0} , {s12*c13,0} , {s13*cdcp,-s13*sdcp} } ,
            { {-s12*c23-c12*s13*s23*cdcp,-c12*s13*s23*sdcp} , {c12*c23-s12*s13*s23*cdcp,-s12*s13*s23*sdcp} , {c13*s23,0} },
            { {s12*s23-c12*s13*c23*cdcp,-c12*s13*c23*sdcp} , {-c12*s23-s12*s13*c23*cdcp,-s12*s13*c23*sdcp} , {c13*c23,0} }
            };

/*
            //print mixing matrix
            std::cout<<""<<std::endl;
            std::cout<<"Mixing matrix U"<<std::endl;
            print_matrix(mixing_matrix);
*/
            //obtain mixin matrix daguer
            std::complex<double> mixing_matrix_daguer[3][3];
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    mixing_matrix_daguer[j][i]=std::conj(mixing_matrix[i][j]);
                };
            };
/*
            //print mixing matrix daguer
            std::cout<<""<<std::endl;
            std::cout<<"Mixing matrix U daguer"<<std::endl;
            print_matrix(mixing_matrix_daguer);
*/


            //prove of U*U_daguer
            std::complex<double> I[3][3];
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    I[i][j]=0;
                    for (int a=0;a<3;a++){
                        I[i][j]=I[i][j]+mixing_matrix[i][a]*mixing_matrix_daguer[a][j];
                    };
                };
            };
/*
            //print I
            std::cout<<""<<std::endl;
            std::cout<<"I matrix"<<std::endl;
            print_matrix(I);
*/

            //mass matrix
            std::complex<double> mass_matrix[3][3];
            mass_matrix[0][0]={m1*m1*c*c*c*c/(2*E),0};
            mass_matrix[1][1]={m2*m2*c*c*c*c/(2*E),0};
            mass_matrix[2][2]={m3*m3*c*c*c*c/(2*E),0};
/*
            //print mass matrix
            std::cout<<""<<std::endl;
            std::cout<<"mass matrix"<<std::endl;
            print_matrix(mass_matrix);
*/

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
/*
            //print vacuum hamiltonian
            std::cout<<""<<std::endl;
            std::cout <<"vacuum hamiltonian"<<std::endl;
            print_matrix(hamiltonian);


            //print initial density matrix
            std::cout<<"···························"<<std::endl;
            std::cout<<"···························"<<std::endl;
            std::cout<<"···························"<<std::endl;
            std::cout <<"initial density matrix"<<std::endl;
            print_matrix(density_matrix);
*/
            
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
/*
            //print k1
            std::cout<<""<<std::endl;
            std::cout <<"k1 matrix"<<std::endl;
            print_matrix(k1);
*/

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
/*
            //print k2
            std::cout<<""<<std::endl;
            std::cout <<"k2 matrix"<<std::endl;
            print_matrix(k2);
*/
  

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
/*
            //print k3
            std::cout<<""<<std::endl;
            std::cout <<"k3 matrix"<<std::endl;
            print_matrix(k3);
*/


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
            
/*             
            //print k4
            std::cout<<""<<std::endl;
            std::cout <<"k4 matrix"<<std::endl;
            print_matrix(k4);
*/            
  

            //found next step density matrix 
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    density_matrix[i][j]=density_matrix[i][j]+(1.0/6.0)*dt*(k1[i][j]+2.0*k2[i][j]+2.0*k3[i][j]+k4[i][j]);
                };
            };  

/*
            //print next step density matrix
            std::cout<<""<<std::endl;
            std::cout <<"final density matrix"<<std::endl;
            print_matrix(density_matrix);
*/

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
};



int main(){

    Density_matrix density;
    std::cout<<real(density.density_matrix[0][0])<<std::endl;    
    for (int i=0;i<100;i++){
        density.evolve_density_matrix(1e-6);
        std::cout<<real(density.density_matrix[0][0])<<std::endl;

    };

    return 0;
};




