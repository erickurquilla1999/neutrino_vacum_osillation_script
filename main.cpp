#include <iostream>
#include <complex>
#include <cmath>

class Density_matrix{
    public:
        //density matrix elements
        double rho_ee;
        double rho_uu;
        double rho_tt;
        double rho_eu;
        double rho_et;
        double rho_ut;

        //initial conditions
        Density_matrix(){
            rho_ee=1;
            rho_et=0;
            rho_tt=0;
            rho_ut=0;
            rho_uu=0;
            rho_eu=0;
        };

        double evolve_density_matrix(){

            //mixing angles radians
            double theta_12=0.23*M_PI;
            double theta_13=0.455*M_PI;
            double theta_23=0.654*M_PI;

            //cp violation angle radians
            double delta_cp=0.633*M_PI;

            //energia IU
            double E=1e-5; 

            //light velocity m/s
            double c=3e-8;

            // sin and cos for U matix
            double s12=sin(theta_12);
            double c12=cos(theta_12);
            double s13=sin(theta_13);
            double c13=cos(theta_13);
            double s23=sin(theta_23);
            double c23=cos(theta_23);
            double sdcp=sin(delta_cp);
            double cdcp=cos(delta_cp);

            //neutrino mases 1=e, 2=u, 3=t IU
            double m1=1;
            double m2=1;
            double m3=1;

            //obtain mixing matrix
            std::complex<double> mixing_matrix[3][3]={
            { {c12*c13,0} , {s12*c13,0} , {s13*cdcp,-s13*sdcp} } ,
            { {-s12*c23-c12*s13*s23*cdcp,-c12*s13*s23*sdcp} , {c12*c23-s12*s13*s23*cdcp,-s12*s13*s23*sdcp} , {c13*s23,0} },
            { {s12*s23-c12*s13*c23*cdcp,-c12*s13*c23*sdcp} , {-c12*s23-s12*s13*c23*cdcp,-s12*s13*c23*sdcp} , {c13*c23,0} }
            };

            //obtain mixin matrix daguer
            std::complex<double> mixing_matrix_daguer[3][3];
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    mixing_matrix_daguer[j][i]=std::conj(mixing_matrix[i][j]);
                };
            };

            //mass matrix
            std::complex<double> mass_matrix[3][3]={{{m1*m1*c*c*c*c/(2*E),0},{0,0},{0,0}},{{0,0},{m2*m2*c*c*c*c/(2*E),0},{0,0}},{{0,0},{0,0},{m3*m3*c*c*c*c/(2*E),0}}};

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

            //print hamiltonian
            std::cout <<"this our hamiltonian"<<std::endl;
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    std::cout << hamiltonian[i][j]<< " ";
                };
                std::cout <<" "<<std::endl;
            };









            return 0.0;

        };
};


int main(){

    Density_matrix density;
    density.evolve_density_matrix();





    return 0;
};








