/*
@author: HuidobroMG


*/

#include <iostream>
#include <math.h>
#include <fstream>
#define M_PI 3.14159265358979323846

// Topological Charge = Baryon Number
int Nb = 1; // The possibilities are: 1, 2, 3, 4, 5, 6, 7, 8 and 4*Nc**3.

// For the cubic B = 4Nc**3 cases
int Nc = cbrt(Nb/4); // Number of B = 4 unit cells
double L = 3.0; // Length of the B = 4 unit cell

// Parameters of the grid
int N = 101; // xmax = (N-1)*dx/2. In the 4Nc**3 cases: N = Nc*(2*L)*1.5/dx + 1
double dx = 0.2;
double dy = 0.2;
double dz = 0.2;
double xmax = (N-1)*dx/2;

// Coefficients of the Lagrangian
double c6 = 0.0; // 2*lambda^2*fpi^2*e^4/hbarc^3
double c0 = 0.2; // 2*mpi^2/(fpi*e)^2

// Step-size for each time evolution
double step_gamma = 2e-5; // if dx = 0.2: gamma = 2e-5; if dx = 0.1: gamma = 2e-6, if dx = 0.4: gamma = 2e-4

// 4 Skyrme fields: 1 Sigma and 3 Pions
int NFields = 4;

// The number of the iteration for the AGD
double iteration = 0.0;

// Initial configuration of the field
void InitPhi(double ****phi0, double *xx, double *yy, double *zz);

// Derivatives functions
double Dx(double ***phi, int i, int j, int k);
double Dy(double ***phi, int i, int j, int k);
double Dz(double ***phi, int i, int j, int k);
double Dxx(double ***phi, int i, int j, int k);
double Dyy(double ***phi, int i, int j, int k);
double Dzz(double ***phi, int i, int j, int k);
double Dxy(double ***phi, int i, int j, int k);
double Dxz(double ***phi, int i, int j, int k);
double Dyz(double ***phi, int i, int j, int k);

// Euler-Lagrange equations
void dEdPhi(double ****phi, double ****psi, double ****phiPrev, double *dB);

// Densities
double BaryonDensity(double ****phi, int i, int j, int k); // The baryon density
double EnergyDensity(double ****phi, int i, int j, int k); // The energy density
double IsospinDensity_11(double ****phi, int i, int j, int k); // The isospin 11 density
double IsospinDensity_22(double ****phi, int i, int j, int k); // The isospin 22 density
double IsospinDensity_33(double ****phi, int i, int j, int k); // The isospin 33 density

// The Baryon Number, Energy and Virial Constraint. It returns the energy
double Baryon_Energy_Virial(double ****phi);
void Isospin(double ****phi); // Isospin inertia tensor

// Derivative of the baryon density
void diffB(double ****phi, int i, int j, int k, double *dB);

// Levi-Civita symbol with 4 indices
double LeviCivita(int i, int j, int k, int l){
    if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)){return 0;}
    else{double eps = static_cast<double>((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l))/12.;
    return eps;}
}

int main(){
    // Create the grid
    double *xx;
    xx = new double [N];
    double *yy;
    yy = new double [N];
    double *zz;
    zz = new double [N];

    for(int i = 0; i < N; i++){
        xx[i] = dx*((double)i - ((double)N-1.0)/2.0);
        yy[i] = dy*((double)i - ((double)N-1.0)/2.0);
        zz[i] = dz*((double)i - ((double)N-1.0)/2.0);
    }

    // Skyrme fields, phi[NFields][N][N][N] and the AGD auxiliar field
    double ****phi;
    double ****psi;
    double ****phiPrev;
    phi = new double ***[NFields]; // phi = [sigma, pi1, pi2, pi3]
    psi = new double ***[NFields];
    phiPrev = new double ***[NFields];
    for(int i=0; i<NFields; i++){
        phi[i] = new double **[N];
        psi[i] = new double **[N];
        phiPrev[i] = new double **[N];
        for(int j=0; j<N; j++){
            phi[i][j] = new double *[N];
            psi[i][j] = new double *[N];
            phiPrev[i][j] = new double *[N];
            for(int k=0; k<N; k++){
                phi[i][j][k] = new double [N];
                psi[i][j][k] = new double [N];
                phiPrev[i][j][k] = new double [N];
            }
        }
    }

    // Initialize the fields
    InitPhi(phi, xx, yy, zz);
    InitPhi(psi, xx, yy, zz);
    InitPhi(phiPrev, xx, yy, zz);
    std::cout << "Ansatz constructed" << std::endl;

    // Define the numerical quantities
    double EPrev = Baryon_Energy_Virial(phi);
    double ENew = 0.0;
    double *dB = new double [3];

    // Start the minimization
    double counter = 0;
    int iter_end = 1000;
    for(int iter = 0; iter < iter_end; iter++){
        std::cout << "iteration = " << iter << " \t" << "velocity = " << iteration << std::endl;
        
        // Update the fields
        dEdPhi(phi, psi, phiPrev, dB);
        iteration += 1.0;
        ENew = Baryon_Energy_Virial(phi);

        // Cut the AGD if the energy increases
        if (EPrev - ENew < 0){iteration = 0.0; counter += 1;}
        EPrev = ENew;
        if (counter == 100){iter = iter_end - 5; counter = 0;}
    }

    // Write the fields in a file
    std::ofstream outfilesigma ("sigma.dat");
    std::ofstream outfilepi1 ("pi1.dat");
    std::ofstream outfilepi2 ("pi2.dat");
    std::ofstream outfilepi3 ("pi3.dat");

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                outfilesigma << phi[0][i][j][k] << " \t";
                outfilepi1 << phi[1][i][j][k] << " \t";
                outfilepi2 << phi[2][i][j][k] << " \t";
                outfilepi3 << phi[3][i][j][k] << " \t";
            }
        }
    }
    outfilesigma.close();
    outfilepi1.close();
    outfilepi2.close();
    outfilepi3.close();

    // Compute the Isospin moment of inertia
    Isospin(phi);
    return 0;
}

// FUNCTIONS
double Baryon_Energy_Virial(double ****phi){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}; // D[coordinate][Field]
    double B = 0.0;
    double E = 0.0;
    double Virial = 0.0;
    double E2 = 0.0;
    double E4 = 0.0;
    double E6 = 0.0;
    double E0 = 0.0;

    for(int i = 2; i<N-2; i++){
        for(int j = 2; j<N-2; j++){
            for(int k = 2; k<N-2; k++){
                for(int a = 0; a<NFields; a++){
                  D[0][a] = Dx(phi[a], i, j, k);
                  D[1][a] = Dy(phi[a], i, j, k);
                  D[2][a] = Dz(phi[a], i, j, k);
                }
                E6 += pow(2.0*M_PI*M_PI*BaryonDensity(phi, i, j, k), 2); // Sextic term
                E0 += 1.0 - phi[0][i][j][k]; // Pion-mass term
                for(int a = 0; a<NFields; a++){
                    E2 += D[0][a]*D[0][a] + D[1][a]*D[1][a] + D[2][a]*D[2][a]; // Quadratic term
                    for(int b = 0; b<NFields; b++){
                        for(int d1 = 0; d1<3; d1++){
                            for(int d2 = 0; d2<3; d2++){
                                E4 += D[d1][a]*D[d1][a]*D[d2][b]*D[d2][b] - D[d1][a]*D[d2][b]*D[d1][b]*D[d2][a]; // Quartic term
                            }
                        }
                        for(int c = 0; c<NFields; c++){
                            for(int d = 0; d<NFields; d++){
                                B += LeviCivita(a,b,c,d)*phi[a][i][j][k]*D[0][b]*D[1][c]*D[2][d];
                            }
                        }
                    }
                }
            }
        }
    }

    B *= -dx*dy*dz/(2.0*M_PI*M_PI);
    E = dx*dy*dz/(24.0*M_PI*M_PI)*(E2 + 2.0*E4 + c6*E6 + c0*E0);
    Virial = dx*dy*dz*abs(1.0 - 2.0*E4/(E2 - 3.0*c6*E6 + 3.0*c0*E0));
    std::cout << "B = " << B << "\t" << "E = " << E << "\t" << "Virial = " << Virial << std::endl;
    return E;
}

void Isospin(double ****phi){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double lambda2_11 = 0;
    double lambda2_22 = 0;
    double lambda2_33 = 0;
    double lambda4_11 = 0;
    double lambda4_22 = 0;
    double lambda4_33 = 0;
    double lambda6_11 = 0;
    double lambda6_22 = 0;
    double lambda6_33 = 0;

    for(int i = 2; i < N-2; i++){
        for(int j = 2; j < N-2; j++){
            for(int k = 2; k < N-2; k++){
                for(int a = 0; a<NFields; a++){
                    D[0][a] = Dx(phi[a],i,j,k);
                    D[1][a] = Dy(phi[a],i,j,k);
                    D[2][a] = Dz(phi[a],i,j,k);
                }
                lambda2_11 += 2*(phi[2][i][j][k]*phi[2][i][j][k] + phi[3][i][j][k]*phi[3][i][j][k]);
                lambda2_22 += 2*(phi[1][i][j][k]*phi[1][i][j][k] + phi[3][i][j][k]*phi[3][i][j][k]);
                lambda2_33 += 2*(phi[1][i][j][k]*phi[1][i][j][k] + phi[2][i][j][k]*phi[2][i][j][k]);
                for (int d1 = 0; d1<3; d1++){
                    lambda4_11 += 8*(D[d1][0]*D[d1][0]*(1-phi[1][i][j][k]*phi[1][i][j][k]) + D[d1][1]*D[d1][1]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[1][i][j][k]*D[d1][0]*D[d1][1]);
                    lambda4_22 += 8*(D[d1][0]*D[d1][0]*(1-phi[2][i][j][k]*phi[2][i][j][k]) + D[d1][2]*D[d1][2]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[2][i][j][k]*D[d1][0]*D[d1][2]);
                    lambda4_33 += 8*(D[d1][0]*D[d1][0]*(1-phi[3][i][j][k]*phi[3][i][j][k]) + D[d1][3]*D[d1][3]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[3][i][j][k]*D[d1][0]*D[d1][3]);
                    for (int d2 = 0; d2<3; d2++){
                        lambda6_11 += (D[d1][0]*D[d2][1] - D[d1][1]*D[d2][0])*(D[d1][0]*D[d2][1] - D[d1][1]*D[d2][0]);
                        lambda6_22 += (D[d1][0]*D[d2][2] - D[d1][2]*D[d2][0])*(D[d1][0]*D[d2][2] - D[d1][2]*D[d2][0]);
                        lambda6_33 += (D[d1][0]*D[d2][3] - D[d1][3]*D[d2][0])*(D[d1][0]*D[d2][3] - D[d1][3]*D[d2][0]);
                    }
                }
            }
        }
    }
    std::cout << "Lambda_11 = " << dx*dy*dz/(24.0*M_PI*M_PI)*(lambda2_11 + lambda4_11 + c6*lambda6_11) << std::endl;
    std::cout << "Lambda_22 = " << dx*dy*dz/(24.0*M_PI*M_PI)*(lambda2_22 + lambda4_22 + c6*lambda6_22) << std::endl;
    std::cout << "Lambda_33 = " << dx*dy*dz/(24.0*M_PI*M_PI)*(lambda2_33 + lambda4_33 + c6*lambda6_33) << std::endl;
}

// Euler-Lagrange equations
void dEdPhi(double ****phi, double ****psi, double ****phiPrev, double *dB){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double DD[3][3][4] = {{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}};
    double dEdphi = 0.0;
    double norma = 1.0;
    double normapsi = 1.0;

    double Vs_0 = 0;
    double Vs_1 = 0;
    double Vs_2 = 0;
    double Lag6_1 = 0;
    double Lag6_2 = 0;

    for(int i = 2; i < N-2; i++){
        for(int j = 2; j < N-2; j++){
            for(int k = 2; k < N-2; k++){
                for(int a = 0; a < NFields; a++){
                    D[0][a] = Dx(phi[a],i,j,k);
                    D[1][a] = Dy(phi[a],i,j,k);
                    D[2][a] = Dz(phi[a],i,j,k);
                    DD[0][0][a] = Dxx(phi[a],i,j,k);
                    DD[1][1][a] = Dyy(phi[a],i,j,k);
                    DD[2][2][a] = Dzz(phi[a],i,j,k);
                    DD[1][0][a] = Dxy(phi[a],i,j,k);
                    DD[0][1][a] = DD[1][0][a];
                    DD[2][0][a] = Dxz(phi[a],i,j,k);
                    DD[0][2][a] = DD[2][0][a];
                    DD[2][1][a] = Dyz(phi[a],i,j,k);
                    DD[1][2][a] = DD[2][1][a];

                    phiPrev[a][i][j][k] = phi[a][i][j][k]; // Save the actual value
                }
                if (c6 != 0 ){
                    Lag6_1 = -2*M_PI*M_PI*BaryonDensity(phi, i, j, k);
                    diffB(phi, i, j, k, dB);
                }
                for(int a = 0; a<NFields; a++){ // There are as many equations as fields we have
                    Lag6_2 = 0;
                    Vs_0 = 0;
                    Vs_1 = 0;
                    Vs_2 = 0;
                    dEdphi = 2*(DD[0][0][a] + DD[1][1][a] + DD[2][2][a]); // Quadratic
                    for(int b = 0; b<NFields; b++){
                        for(int d1 = 0; d1<3; d1++){
                            for(int d2 = 0; d2<3; d2++){
                                dEdphi += 8*(2*DD[d1][d2][b]*D[d1][a]*D[d2][b] + DD[d1][d1][a]*D[d2][b]*D[d2][b] - DD[d1][d2][a]*D[d2][b]*D[d1][b] - DD[d1][d2][b]*D[d2][a]*D[d1][b] - DD[d1][d1][b]*D[d2][a]*D[d2][b]); // Quartic
                            }
                        }
                        if (c6 != 0){
                           for(int c = 0; c<NFields; c++){
                                for(int d = 0; d<NFields; d++){
                                    Lag6_2 += LeviCivita(a,b,c,d)*D[0][b]*D[1][c]*D[2][d];
                                    Vs_0 -= LeviCivita(a,b,c,d)*phi[b][i][j][k]*D[1][c]*D[2][d];
                                    Vs_1 += LeviCivita(a,b,c,d)*phi[b][i][j][k]*D[0][c]*D[2][d];
                                    Vs_2 -= LeviCivita(a,b,c,d)*phi[b][i][j][k]*D[0][c]*D[1][d];
                                }
                            }
                        }
                    }
                    if (c6 != 0){dEdphi += -8*c6*Lag6_1*Lag6_2 + 2*c6*(Vs_0*dB[0] + Vs_1*dB[1] + Vs_2*dB[2]);} // Sextic
                    phi[a][i][j][k] = psi[a][i][j][k] + step_gamma*dEdphi + step_gamma*(a == 0)*c0; // Accelerated Gradient Descent
                }

                // We have to impose the normalization condition of the fields
                norma = sqrt(phi[0][i][j][k]*phi[0][i][j][k] + phi[1][i][j][k]*phi[1][i][j][k] + phi[2][i][j][k]*phi[2][i][j][k] + phi[3][i][j][k]*phi[3][i][j][k]);
                for(int a = 0; a < NFields; a++){
                  phi[a][i][j][k] /= norma;
                  psi[a][i][j][k] = phi[a][i][j][k] + iteration/(iteration + 3.0)*(phi[a][i][j][k] - phiPrev[a][i][j][k]);
                }

                normapsi = sqrt(psi[0][i][j][k]*psi[0][i][j][k] + psi[1][i][j][k]*psi[1][i][j][k] + psi[2][i][j][k]*psi[2][i][j][k] + psi[3][i][j][k]*psi[3][i][j][k]);
                for(int a = 0; a < NFields; a++){
                  psi[a][i][j][k] /= normapsi;
                }
            }
        }
    }
    /*
    // Vacuum Boundary Conditions
    for(int a = 0; a < NFields; a++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < N; k++){
                phi[a][0][j][k] = 0.0;
                phi[a][1][j][k] = 0.0;
                phi[a][N-2][j][k] = 0.0;
                phi[a][N-1][j][k] = 0.0;

                phi[a][j][0][k] = 0.0;
                phi[a][j][1][k] = 0.0;
                phi[a][j][N-2][k] = 0.0;
                phi[a][j][N-1][k] = 0.0;

                phi[a][j][k][0] = 0.0;
                phi[a][j][k][1] = 0.0;
                phi[a][j][k][N-2] = 0.0;
                phi[a][j][k][N-1] = 0.0;

                psi[a][0][j][k] = 0.0;
                psi[a][1][j][k] = 0.0;
                psi[a][N-2][j][k] = 0.0;
                psi[a][N-1][j][k] = 0.0;

                psi[a][j][0][k] = 0.0;
                psi[a][j][1][k] = 0.0;
                psi[a][j][N-2][k] = 0.0;
                psi[a][j][N-1][k] = 0.0;

                psi[a][j][k][0] = 0.0;
                psi[a][j][k][1] = 0.0;
                psi[a][j][k][N-2] = 0.0;
                psi[a][j][k][N-1] = 0.0;
            }
        }
    }
    */
}

// Derivative of the Baryon density (Auxiliar function for faster performance)
void diffB(double ****phi, int i, int j, int k, double *dB){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double DD[3][3][4] = {{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}};

    dB[0] = 0;
    dB[1] = 0;
    dB[2] = 0;
    for(int a = 0; a<NFields; a++){
      D[0][a] = Dx(phi[a], i, j, k);
      D[1][a] = Dy(phi[a], i, j, k);
      D[2][a] = Dz(phi[a], i, j, k);
      DD[0][0][a] = Dxx(phi[a],i,j,k);
      DD[1][1][a] = Dyy(phi[a],i,j,k);
      DD[2][2][a] = Dzz(phi[a],i,j,k);
      DD[1][0][a] = Dxy(phi[a],i,j,k);
      DD[0][1][a] = DD[1][0][a];
      DD[2][0][a] = Dxz(phi[a],i,j,k);
      DD[0][2][a] = DD[2][0][a];
      DD[2][1][a] = Dyz(phi[a],i,j,k);
      DD[1][2][a] = DD[2][1][a];
    }
    for(int a = 0; a<NFields; a++){
        for(int b = 0; b<NFields; b++){
            for(int c = 0; c<NFields; c++){
                for(int d = 0; d<NFields; d++){
                    dB[0] += LeviCivita(a,b,c,d)*phi[a][i][j][k]*(DD[0][0][b]*D[1][c]*D[2][d] + DD[0][1][c]*D[0][b]*D[2][d] + DD[0][2][d]*D[0][b]*D[1][c]);
                    dB[1] += LeviCivita(a,b,c,d)*phi[a][i][j][k]*(DD[1][0][b]*D[1][c]*D[2][d] + DD[1][1][c]*D[0][b]*D[2][d] + DD[1][2][d]*D[0][b]*D[1][c]);
                    dB[2] += LeviCivita(a,b,c,d)*phi[a][i][j][k]*(DD[2][0][b]*D[1][c]*D[2][d] + DD[2][1][c]*D[0][b]*D[2][d] + DD[2][2][d]*D[0][b]*D[1][c]);
                }
            }
        }
    }
}

// Rational Map Ansatz
void InitPhi(double ****phi0, double *xx, double *yy, double *zz){
    // Rational map ansatz
    double zRe = 0;
    double zIm = 0;
    double RRe = 0;
    double RIm = 0;
    double Rmod2 = 0;
    double Denom = 1;

    double sign = 1;
    double sig_p = 0;
    double pi1_p = 0;
    double pi2_p = 0;
    double pi3_p = 0;
    double sig = 0;
    double pi1 = 0;
    double pi2 = 0;
    double pi3 = 0;
    double t = 0;
    double qx = 0;
    double qy = 0;
    double qz = 0;
    double qr = 0;
    double a_out = 0;
    double b_out = 0;
    double c_out = 0;

    double r = 0;
    double rho = 0;
    double f = 0;
    double theta = 0;
    double pangle = 0;

    if (Nb <= 8){
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                for(int k = 0; k<N; k++){
                    r = sqrt(xx[i]*xx[i] + yy[j]*yy[j] + zz[k]*zz[k]);
                    rho = sqrt(xx[i]*xx[i] + yy[j]*yy[j]);
                    if(c0 != 0){f = M_PI*exp(-r/(xmax/7.0));}
                    else{f = M_PI/(1.0 + r*r/(xmax/5.0));}
                    theta = atan2(rho, zz[k]);
                    pangle = atan2(yy[j], xx[i]);

                    zRe = tan(theta/2.0)*cos(pangle);
                    zIm = tan(theta/2.0)*sin(pangle);

                    if (Nb == 1){
                        RRe = zRe;
                        RIm = zIm;
                    }
                    else if (Nb == 2){
                        RRe = zRe*zRe - zIm*zIm;
                        RIm = 2*zRe*zIm;
                    }
                    else if (Nb == 3){
                        Denom = pow(zIm, 6) + 3*pow(zIm, 4)*zRe*zRe - 4*sqrt(3)*pow(zIm, 3)*zRe + 3*zIm*zIm*pow(zRe, 4) + 3*zIm*zIm - 4*sqrt(3)*zIm*pow(zRe, 3) + pow(zRe, 6) + 3*pow(zRe, 2);
                        if(Denom == 0){Denom = 1e-15;}
                        RRe = (sqrt(3)*pow(zIm, 5) + 2*sqrt(3)*pow(zIm, 3)*pow(zRe, 2) + sqrt(3)*zIm*pow(zRe, 4) - sqrt(3)*zIm - 4*pow(zRe, 3))/Denom;
                        RIm = (sqrt(3)*pow(zIm, 4)*zRe - 4*pow(zIm, 3) + 2*sqrt(3)*pow(zIm, 2)*pow(zRe, 3) + sqrt(3)*pow(zRe, 5) - sqrt(3)*zRe)/Denom;
                    }
                    else if (Nb == 4){
                        Denom = ((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm + 4*sqrt(3.0)*zRe*zIm + 1.0)*((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm + 4*sqrt(3.0)*zRe*zIm + 1.0) + (4*zRe*zIm*(zRe*zRe - zIm*zIm) - 2*sqrt(3.0)*(zRe*zRe - zIm*zIm))*(4*zRe*zIm*(zRe*zRe - zIm*zIm) - 2*sqrt(3.0)*(zRe*zRe - zIm*zIm));
                        RRe = (((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm - 4*sqrt(3.0)*zRe*zIm + 1.0)*((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm + 4*sqrt(3.0)*zRe*zIm + 1.0) + (4*zRe*zIm*(zRe*zRe - zIm*zIm) + 2*sqrt(3.0)*(zRe*zRe - zIm*zIm))*(4*zRe*zIm*(zRe*zRe - zIm*zIm) - 2*sqrt(3.0)*(zRe*zRe - zIm*zIm)))/Denom;
                        RIm = (-((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm - 4*sqrt(3.0)*zRe*zIm + 1.0)*(4*zRe*zIm*(zRe*zRe - zIm*zIm) - 2*sqrt(3.0)*(zRe*zRe - zIm*zIm)) + (4*zRe*zIm*(zRe*zRe - zIm*zIm) + 2*sqrt(3.0)*(zRe*zRe - zIm*zIm))*((zRe*zRe - zIm*zIm)*(zRe*zRe - zIm*zIm) - 4*zRe*zRe*zIm*zIm + 4*sqrt(3.0)*zRe*zIm + 1.0))/Denom;
                    }
                    else if (Nb == 5){
                        Denom = 9.42*pow(zIm,8) + 37.7*pow(zIm,6)*pow(zRe,2) + 24.19*pow(zIm,6) + 56.55*pow(zIm,4)*pow(zRe,4) + 24.19*pow(zIm,4)*pow(zRe,2) + 21.66*pow(zIm,4) + 37.7*pow(zIm,2)*pow(zRe,6) - 24.19*pow(zIm,2)*pow(zRe,4) - 5.79*pow(zIm,2)*pow(zRe,2) + 7.88*pow(zIm,2) + 9.42*pow(zRe,8) - 24.19*pow(zRe,6) + 21.66*pow(zRe,4) - 7.88*pow(zRe,2) + 1.0;
                        RRe = (3.07*pow(zIm,8)*zRe + 12.28*pow(zIm,6)*pow(zRe,3) + 23.92*pow(zIm,6)*zRe + 18.42*pow(zIm,4)*pow(zRe,5) + 55.99*pow(zIm,4)*pow(zRe,3) - 38.8*pow(zIm,4)*zRe + 12.28*pow(zIm,2)*pow(zRe,7) + 40.23*pow(zIm,2)*pow(zRe,5) - 59.9*pow(zIm,2)*pow(zRe,3) - 23.92*pow(zIm,2)*zRe + 3.07*pow(zRe,9) + 8.16*pow(zRe,7) - 5.1*pow(zRe,5) - 8.16*pow(zRe,3) + 3.07*zRe)/Denom;
                        RIm = (3.07*pow(zIm,9) + 12.28*pow(zIm,7)*pow(zRe,2) - 8.16*pow(zIm,7) + 18.42*pow(zIm,5)*pow(zRe,4) - 40.23*pow(zIm,5)*pow(zRe,2) - 5.1*pow(zIm,5) + 12.28*pow(zIm,3)*pow(zRe,6) - 55.99*pow(zIm,3)*pow(zRe,4) - 59.9*pow(zIm,3)*pow(zRe,2) + 8.16*pow(zIm,3) + 3.07*zIm*pow(zRe,8) - 23.92*zIm*pow(zRe,6) - 38.8*zIm*pow(zRe,4) + 23.92*zIm*pow(zRe,2) + 3.07*zIm)/Denom;
                    }
                    else if (Nb == 6){
                        Denom = 0.03*pow(zIm,12) + 0.15*pow(zIm,10)*pow(zRe,2) + 0.38*pow(zIm,8)*pow(zRe,4) + 1.28*pow(zIm,7)*zRe + 0.51*pow(zIm,6)*pow(zRe,6) + 1.28*pow(zIm,5)*pow(zRe,3) + 0.38*pow(zIm,4)*pow(zRe,8) + pow(zIm,4) - 1.28*pow(zIm,3)*pow(zRe,5) + 0.15*pow(zIm,2)*pow(zRe,10) + 2.0*pow(zIm,2)*pow(zRe,2) - 1.28*zIm*pow(zRe,7) + 0.03*pow(zRe,12) + pow(zRe,4);
                        if(Denom == 0){Denom = 1e-15;}
                        RRe = (-0.32*pow(zIm,9)*zRe - 1.28*pow(zIm,7)*pow(zRe,3) - 1.03*pow(zIm,6) - 1.92*pow(zIm,5)*pow(zRe,5) - 0.62*pow(zIm,4)*pow(zRe,2) - 1.28*pow(zIm,3)*pow(zRe,7) + 0.62*pow(zIm,2)*pow(zRe,4) - 0.32*zIm*pow(zRe,9) + 0.32*zIm*zRe + 1.03*pow(zRe,6))/Denom;
                        RIm = (0.16*pow(zIm,10) + 0.48*pow(zIm,8)*pow(zRe,2) + 0.32*pow(zIm,6)*pow(zRe,4) + 1.85*pow(zIm,5)*zRe - 0.32*pow(zIm,4)*pow(zRe,6) + 4.51*pow(zIm,3)*pow(zRe,3) - 0.48*pow(zIm,2)*pow(zRe,8) - 0.16*pow(zIm,2) + 1.85*zIm*pow(zRe,5) - 0.16*pow(zRe,10) + 0.16*pow(zRe,2))/Denom;
                    }
                    else if (Nb == 7){
                        Denom = 0.02*pow(zIm,14) + 0.14*pow(zIm,12)*zRe*zRe + 0.43*pow(zIm,10)*pow(zRe,4) + 0.71*pow(zIm,8)*pow(zRe,6) - 1.43*pow(zIm,8)*zRe + 0.71*pow(zIm,6)*pow(zRe,8) + 0.43*pow(zIm,4)*pow(zRe,10) + 4*pow(zIm,4)*pow(zRe,5) + pow(zIm,4) + 0.14*pow(zIm,2)*pow(zRe,12) + 2.29*pow(zIm,2)*pow(zRe,7) + 2*zIm*zIm*zRe*zRe + 0.02*pow(zRe,14) - 0.29*pow(zRe,9) + pow(zRe,4);
                        if(Denom == 0){Denom = 1e-15;}
                        RRe = (0.14*pow(zIm,12) + 0.57*pow(zIm,10)*zRe*zRe + 0.71*pow(zIm,8)*pow(zRe,4) - 2.86*pow(zIm,6)*zRe - 0.71*pow(zIm,4)*pow(zRe,8) - 5.71*pow(zIm,4)*pow(zRe,3) - 0.57*pow(zIm,2)*pow(zRe,10) - 0.57*pow(zIm,2)*pow(zRe,5) - 0.14*zIm*zIm - 0.14*pow(zRe,12) + 0.98*pow(zRe,7) + 0.14*zRe*zRe)/Denom;
                        RIm = (0.29*pow(zIm,11)*zRe + 1.43*pow(zIm,9)*pow(zRe,3) + 2.86*pow(zIm,7)*pow(zRe,5) - 1.02*pow(zIm,7) + 2.86*pow(zIm,5)*pow(zRe,7) + 1.43*pow(zIm,5)*zRe*zRe + 1.43*pow(zIm,3)*pow(zRe,9) + 4.29*pow(zIm,3)*pow(zRe,4) + 0.29*zIm*pow(zRe,11) + 3.14*zIm*pow(zRe,6) - 0.29*zIm*zRe)/Denom;
                    }
                    else if (Nb == 8){
                        Denom = (pow(zRe,6) - 15*pow(zRe,4)*zIm*zIm + 15*zRe*zRe*pow(zIm,4) - pow(zIm,6) - 0.14)*(pow(zRe,6) - 15*pow(zRe,4)*zIm*zIm + 15*zRe*zRe*pow(zIm,4) - pow(zIm,6) - 0.14) + (6*zIm*pow(zRe,5) - 20*pow(zRe*zIm,3) + 6*zRe*pow(zIm,5))*(6*zIm*pow(zRe,5) - 20*pow(zRe*zIm,3) + 6*zRe*pow(zIm,5));
                        RRe = ((0.14*pow(zRe,8) - 3.92*zIm*zIm*pow(zRe,6) + 9.8*pow(zRe*zIm,4) - 3.92*zRe*zRe*pow(zIm,6) + zRe*zRe - zIm*zIm + 0.14*pow(zIm,8))*(pow(zRe,6) - 15*pow(zRe,4)*zIm*zIm + 15*zRe*zRe*pow(zIm,4) - pow(zIm,6) - 0.14) + (1.12*zIm*pow(zRe,7) - 7.84*pow(zIm,3)*pow(zRe,5) + 7.84*pow(zRe,3)*pow(zIm,5) - 1.12*zRe*pow(zIm,7) + 2*zRe*zIm)*(6*zIm*pow(zRe,5) - 20*pow(zRe*zIm,3) + 6*zRe*pow(zIm,5)))/Denom;
                        RIm = (-(0.14*pow(zRe,8) - 3.92*zIm*zIm*pow(zRe,6) + 9.8*pow(zRe*zIm,4) - 3.92*zRe*zRe*pow(zIm,6) + zRe*zRe - zIm*zIm + 0.14*pow(zIm,8))*(6*zIm*pow(zRe,5) - 20*pow(zRe*zIm,3) + 6*zRe*pow(zIm,5)) + (1.12*zIm*pow(zRe,7) - 7.84*pow(zIm,3)*pow(zRe,5) + 7.84*pow(zRe,3)*pow(zIm,5) - 1.12*zRe*pow(zIm,7) + 2*zRe*zIm)*(pow(zRe,6) - 15*pow(zRe,4)*zIm*zIm + 15*zRe*zRe*pow(zIm,4) - pow(zIm,6) - 0.14))/Denom;
                    }

                    Rmod2 = RRe*RRe + RIm*RIm;
                    phi0[0][i][j][k] = cos(f);
                    phi0[1][i][j][k] = sin(f)*2*RRe/(1.0 + Rmod2);
                    phi0[2][i][j][k] = sin(f)*2*RIm/(1.0 + Rmod2);
                    phi0[3][i][j][k] = sin(f)*(1.0 - Rmod2)/(1.0 + Rmod2);
                }
            }
        }
    }
    else{
        if(Nc % 2 != 0.0){sign = -1;} // Even or Odd
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                for(int k = 0; k<N; k++){
                    r = sqrt(xx[i]*xx[i] + yy[j]*yy[j] + zz[k]*zz[k]);
                    if (abs(xx[i]) <= Nc*L && abs(yy[j]) <= Nc*L && abs(zz[k]) <= Nc*L){ // Inside the Crystal
                        // Construct the fields with the Castillejo ansatz
                        sig = sin(M_PI*xx[i]/L)*sin(M_PI*yy[j]/L)*sin(M_PI*zz[k]/L);
                        pi1 = cos(M_PI*xx[i]/L)*sqrt(1.0 - 0.5*cos(M_PI*yy[j]/L)*cos(M_PI*yy[j]/L) - 0.5*cos(M_PI*zz[k]/L)*cos(M_PI*zz[k]/L) + 1.0/3.0*cos(M_PI*yy[j]/L)*cos(M_PI*yy[j]/L)*cos(M_PI*zz[k]/L)*cos(M_PI*zz[k]/L));
                        pi2 = cos(M_PI*yy[j]/L)*sqrt(1.0 - 0.5*cos(M_PI*zz[k]/L)*cos(M_PI*zz[k]/L) - 0.5*cos(M_PI*xx[i]/L)*cos(M_PI*xx[i]/L) + 1.0/3.0*cos(M_PI*zz[k]/L)*cos(M_PI*zz[k]/L)*cos(M_PI*xx[i]/L)*cos(M_PI*xx[i]/L));
                        pi3 = cos(M_PI*zz[k]/L)*sqrt(1.0 - 0.5*cos(M_PI*xx[i]/L)*cos(M_PI*xx[i]/L) - 0.5*cos(M_PI*yy[j]/L)*cos(M_PI*yy[j]/L) + 1.0/3.0*cos(M_PI*xx[i]/L)*cos(M_PI*xx[i]/L)*cos(M_PI*yy[j]/L)*cos(M_PI*yy[j]/L));

                        // Perform the O(4) isorotation to link with the vacuum
                        phi0[0][i][j][k] = sign/sqrt(3)*(pi1 + pi2 + pi3);
                        phi0[1][i][j][k] = -sign/sqrt(3)*sig + 1.0/3.0*(2*pi1 - pi2 - pi3);
                        phi0[2][i][j][k] = -sign/sqrt(3)*sig + 1.0/3.0*(2*pi2 - pi3 - pi1);
                        phi0[3][i][j][k] = -sign/sqrt(3)*sig + 1.0/3.0*(2*pi3 - pi1 - pi2);
                    }
                    else{ // Outside the Crystal
                        if (abs(xx[i]) >= abs(yy[j]) && abs(xx[i]) >= abs(zz[k])){
                            if (xx[i]<0){t = -Nc*L/xx[i];}
                            else{t = Nc*L/xx[i];}
                        }
                        else if (abs(yy[j])>=abs(xx[i]) && abs(yy[j])>=abs(zz[k])){
                            if (yy[j]<0){t = -Nc*L/yy[j];}
                            else{t = Nc*L/yy[j];}
                        }
                        else if (abs(zz[k])>=abs(xx[i]) && abs(zz[k])>=abs(yy[j])){
                            if (zz[k]<0){t = -Nc*L/zz[k];}
                            else{t = Nc*L/zz[k];}
                        }
                        // Project the coordinates in the nearest surface
                        qx = t*xx[i];
                        qy = t*yy[j];
                        qz = t*zz[k];
                        qr = sqrt(qx*qx + qy*qy + qz*qz);

                        // Construct the fields in that surface
                        sig_p = sin(M_PI*qx/L)*sin(M_PI*qy/L)*sin(M_PI*qz/L);
                        pi1_p = cos(M_PI*qx/L)*sqrt(1.0 - 0.5*cos(M_PI*qy/L)*cos(M_PI*qy/L) - 0.5*cos(M_PI*qz/L)*cos(M_PI*qz/L) + 1.0/3.0*cos(M_PI*qy/L)*cos(M_PI*qy/L)*cos(M_PI*qz/L)*cos(M_PI*qz/L));
                        pi2_p = cos(M_PI*qy/L)*sqrt(1.0 - 0.5*cos(M_PI*qz/L)*cos(M_PI*qz/L) - 0.5*cos(M_PI*qx/L)*cos(M_PI*qx/L) + 1.0/3.0*cos(M_PI*qz/L)*cos(M_PI*qz/L)*cos(M_PI*qx/L)*cos(M_PI*qx/L));
                        pi3_p = cos(M_PI*qz/L)*sqrt(1.0 - 0.5*cos(M_PI*qx/L)*cos(M_PI*qx/L) - 0.5*cos(M_PI*qy/L)*cos(M_PI*qy/L) + 1.0/3.0*cos(M_PI*qx/L)*cos(M_PI*qx/L)*cos(M_PI*qy/L)*cos(M_PI*qy/L));

                        // Perform the O(4) isorotation to link with the vacuum
                        sig = sign/sqrt(3.000001)*(pi1_p + pi2_p + pi3_p);
                        pi1 = -sign/sqrt(3.000001)*sig_p + 1.0/3.000001*(2*pi1_p - pi2_p - pi3_p);
                        pi2 = -sign/sqrt(3.000001)*sig_p + 1.0/3.000001*(2*pi2_p - pi3_p - pi1_p);
                        pi3 = -sign/sqrt(3.000001)*sig_p + 1.0/3.000001*(2*pi3_p - pi1_p - pi2_p);

                        a_out = acos(sig);
                        b_out = acos(pi1/sin(a_out));
                        if (pi3 >= 0){c_out = acos(pi2/(sin(a_out)*sin(b_out)));}
                        else{c_out = 2*M_PI - acos(pi2/(sin(a_out)*sin(b_out)));}

                        // Construct the fields interpolating the vacuum
                        phi0[0][i][j][k] = cos(exp(-r/qr)*a_out);
                        phi0[1][i][j][k] = sin(exp(-r/qr)*a_out)*cos(b_out);
                        phi0[2][i][j][k] = sin(exp(-r/qr)*a_out)*sin(b_out)*cos(c_out);
                        phi0[3][i][j][k] = sin(exp(-r/qr)*a_out)*sin(b_out)*sin(c_out);
                    }
                }
            }
        }
    }
}

// Define the derivatives
double Dx(double ***phi, int i, int j, int k){
    return (-phi[i+2][j][k] + 8.0*phi[i+1][j][k] - 8.0*phi[i-1][j][k] + phi[i-2][j][k])/(12.0*dx);
}

double Dy(double ***phi, int i, int j, int k){
    return (-phi[i][j+2][k] + 8.0*phi[i][j+1][k] - 8.0*phi[i][j-1][k] + phi[i][j-2][k])/(12.0*dy);
}

double Dz(double ***phi, int i, int j, int k){
    return (-phi[i][j][k+2] + 8.0*phi[i][j][k+1] - 8.0*phi[i][j][k-1] + phi[i][j][k-2])/(12.0*dz);
}

double Dxx(double ***phi, int i, int j, int k){
    return (-phi[i+2][j][k] + 16.0*phi[i+1][j][k] - 30.0*phi[i][j][k] + 16.0*phi[i-1][j][k] - phi[i-2][j][k])/(12.0*dx*dx);
}

double Dyy(double ***phi, int i, int j, int k){
    return (-phi[i][j+2][k] + 16.0*phi[i][j+1][k] - 30.0*phi[i][j][k] + 16.0*phi[i][j-1][k] - phi[i][j-2][k])/(12.0*dy*dy);
}

double Dzz(double ***phi, int i, int j, int k){
    return (-phi[i][j][k+2] + 16.0*phi[i][j][k+1] - 30.0*phi[i][j][k] + 16.0*phi[i][j][k-1] - phi[i][j][k-2])/(12.0*dz*dz);
}

double Dxy(double ***phi, int i, int j, int k){
    return (phi[i+2][j][k] + phi[i-2][j][k] + phi[i][j-2][k] + phi[i][j+2][k] - phi[i+2][j+2][k] - phi[i-2][j-2][k] - 16.0*phi[i+1][j][k] + 30.0*phi[i][j][k] - 16.0*phi[i-1][j][k] - 16.0*phi[i][j+1][k] - 16.0*phi[i][j-1][k] + 16.0*phi[i+1][j+1][k] + 16.0*phi[i-1][j-1][k])/(24.0*dx*dy);
}

double Dxz(double ***phi, int i, int j, int k){
    return (phi[i+2][j][k] + phi[i-2][j][k] + phi[i][j][k-2] + phi[i][j][k+2] - phi[i+2][j][k+2] - phi[i-2][j][k-2] - 16.0*phi[i+1][j][k] + 30.0*phi[i][j][k] - 16.0*phi[i-1][j][k] - 16.0*phi[i][j][k+1] - 16.0*phi[i][j][k-1] + 16.0*phi[i+1][j][k+1] + 16.0*phi[i-1][j][k-1])/(24.0*dx*dz);
}

double Dyz(double ***phi, int i, int j, int k){
    return (phi[i][j+2][k] + phi[i][j-2][k] + phi[i][j][k-2] + phi[i][j][k+2] - phi[i][j+2][k+2] - phi[i][j-2][k-2] - 16.0*phi[i][j+1][k] + 30.0*phi[i][j][k] - 16.0*phi[i][j-1][k] - 16.0*phi[i][j][k+1] - 16.0*phi[i][j][k-1] + 16.0*phi[i][j+1][k+1] + 16.0*phi[i][j-1][k-1])/(24.0*dy*dz);
}

double BaryonDensity(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double Bd = 0;

    for(int a = 0; a<NFields; a++){
      D[0][a] = Dx(phi[a], i, j, k);
      D[1][a] = Dy(phi[a], i, j, k);
      D[2][a] = Dz(phi[a], i, j, k);
    }
    for(int a = 0; a<NFields; a++){
        for(int b = 0; b<NFields; b++){
            for(int c = 0; c<NFields; c++){
                for(int d = 0; d<NFields; d++){
                    Bd += LeviCivita(a,b,c,d)*phi[a][i][j][k]*D[0][b]*D[1][c]*D[2][d];
                }
            }
        }
    }

    return -Bd/(2*M_PI*M_PI);
}

double EnergyDensity(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double E2 = 0;
    double E4 = 0;
    double E6 = 0;
    double E0 = 0;

    E6 = pow(2*M_PI*M_PI*BaryonDensity(phi, i, j, k), 2); // Sextic term
    E0 = 1 - phi[0][i][j][k]; // Pion-mass term
    for(int a = 0; a<NFields; a++){
        D[0][a] = Dx(phi[a],i,j,k);
        D[1][a] = Dy(phi[a],i,j,k);
        D[2][a] = Dz(phi[a],i,j,k);
    }
    for(int a = 0; a<NFields; a++){
        E2 += D[0][a]*D[0][a] + D[1][a]*D[1][a] + D[2][a]*D[2][a]; // Quadratic term
        for(int b = 0; b<NFields; b++){
            for(int d1 = 0; d1<3; d1++){
                for(int d2 = 0; d2<3; d2++){
                    E4 += D[d1][a]*D[d1][a]*D[d2][b]*D[d2][b] - D[d1][a]*D[d2][b]*D[d1][b]*D[d2][a]; // Quartic term
                }
            }
        }
    }
    return 1.0/(24.0*M_PI*M_PI)*(E2 + 2.0*E4 + c6*E6 + c0*E0);;
}

double IsospinDensity_11(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double lambda2_11 = 0;
    double lambda4_11 = 0;
    double lambda6_11 = 0;

    for(int i = 2; i < N-2; i++){
        for(int j = 2; j < N-2; j++){
            for(int k = 2; k < N-2; k++){
                for(int a = 0; a<NFields; a++){
                    D[0][a] = Dx(phi[a],i,j,k);
                    D[1][a] = Dy(phi[a],i,j,k);
                    D[2][a] = Dz(phi[a],i,j,k);
                }
                lambda2_11 += 2*(phi[2][i][j][k]*phi[2][i][j][k] + phi[3][i][j][k]*phi[3][i][j][k]);
                for (int d1 = 0; d1<3; d1++){
                    lambda4_11 += 8*(D[d1][0]*D[d1][0]*(1-phi[1][i][j][k]*phi[1][i][j][k]) + D[d1][1]*D[d1][1]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[1][i][j][k]*D[d1][0]*D[d1][1]);
                    for (int d2 = 0; d2<3; d2++){
                        lambda6_11 += (D[d1][0]*D[d2][1] - D[d1][1]*D[d2][0])*(D[d1][0]*D[d2][1] - D[d1][1]*D[d2][0]);
                    }
                }
            }
        }
    }
    return 1.0/(24.0*M_PI*M_PI)*(lambda2_11 + lambda4_11 + c6*lambda6_11);
}

double IsospinDensity_22(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double lambda2_22 = 0;
    double lambda4_22 = 0;
    double lambda6_22 = 0;

    for(int i = 2; i < N-2; i++){
        for(int j = 2; j < N-2; j++){
            for(int k = 2; k < N-2; k++){
                for(int a = 0; a<NFields; a++){
                    D[0][a] = Dx(phi[a],i,j,k);
                    D[1][a] = Dy(phi[a],i,j,k);
                    D[2][a] = Dz(phi[a],i,j,k);
                }
                lambda2_22 += 2*(phi[1][i][j][k]*phi[1][i][j][k] + phi[3][i][j][k]*phi[3][i][j][k]);
                for (int d1 = 0; d1<3; d1++){
                    lambda4_22 += 8*(D[d1][0]*D[d1][0]*(1-phi[2][i][j][k]*phi[2][i][j][k]) + D[d1][2]*D[d1][2]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[2][i][j][k]*D[d1][0]*D[d1][2]);
                    for (int d2 = 0; d2<3; d2++){
                        lambda6_22 += (D[d1][0]*D[d2][2] - D[d1][2]*D[d2][0])*(D[d1][0]*D[d2][2] - D[d1][2]*D[d2][0]);
                    }
                }
            }
        }
    }
    return 1.0/(24.0*M_PI*M_PI)*(lambda2_22 + lambda4_22 + c6*lambda6_22);
}

double IsospinDensity_33(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double lambda2_33 = 0;
    double lambda4_33 = 0;
    double lambda6_33 = 0;

    for(int i = 2; i < N-2; i++){
        for(int j = 2; j < N-2; j++){
            for(int k = 2; k < N-2; k++){
                for(int a = 0; a<NFields; a++){
                    D[0][a] = Dx(phi[a],i,j,k);
                    D[1][a] = Dy(phi[a],i,j,k);
                    D[2][a] = Dz(phi[a],i,j,k);
                }
                lambda2_33 += 2*(phi[1][i][j][k]*phi[1][i][j][k] + phi[2][i][j][k]*phi[2][i][j][k]);
                for (int d1 = 0; d1<3; d1++){
                    lambda4_33 += 8*(D[d1][0]*D[d1][0]*(1-phi[3][i][j][k]*phi[3][i][j][k]) + D[d1][3]*D[d1][3]*(1-phi[0][i][j][k]*phi[0][i][j][k]) + 2*phi[0][i][j][k]*phi[3][i][j][k]*D[d1][0]*D[d1][3]);
                    for (int d2 = 0; d2<3; d2++){
                        lambda6_33 += (D[d1][0]*D[d2][3] - D[d1][3]*D[d2][0])*(D[d1][0]*D[d2][3] - D[d1][3]*D[d2][0]);
                    }
                }
            }
        }
    }
    return 1.0/(24.0*M_PI*M_PI)*(lambda2_33 + lambda4_33 + c6*lambda6_33);
}

