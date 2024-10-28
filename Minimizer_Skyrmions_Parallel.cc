/*
@author: HuidobroMG

This code is the paralellized version of Minimizer_Skyrmions.cc using the openmp library of C++.
It must be compiled using,
g++ -fopenmp output.exe Minimizer_Skyrmions_Parallel.cc
*/

#include <iostream>
#include <math.h>
#include <fstream>
#include <omp.h>
#define M_PI 3.14159265358979323846

// Number of threads for parallel jobs
#define numThreads 14 // Use omp_get_max_threads() for max threads

// Parameters of the grid
int N = 101; // xmax = (N-1)*dx/2. In the 4Nc**3 cases: N = Nc*2*L*1.5/dx + 1
double dx = 0.2;
double dy = 0.2;
double dz = 0.2;
double xmax = (N-1)*dx/2;

// Coefficients of the Lagrangian
double c6 = 0.0; // 2*lambda^2*fpi^2*e^4/hbarc^3
double c0 = 0.2; // 2*mpi^2/(fpi*e)^2

// Step-size for each time evolution
double step_gamma = 2e-5; // if dx = 0.2, gamma = 2e-5; if dx = 0.1, gamma = 2e-6, if dx = 0.4, gamma = 2e-4

// 4 Skyrme fields: 1 Sigma and 3 Pions
int NFields = 4;

// The number of the iteration for the AGD
double iteration = 0.0;

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

// The integrals
double Baryon(double ****phi);
double Energy(double ****phi);
double Virial(double ****phi);
void Isospin(double ****phi);

// Vector derivative of the baryon density
void diffB(double ****phi, int i, int j, int k, double *dB);

// Levi-Civita symbol with 4 indices
double LeviCivita(int i, int j, int k, int l){
    if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)){return 0;}
    else{double eps = static_cast<double>((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l))/12.;
    return eps;}
}

int main(){
    // Set number of threads for parallel regions
    omp_set_dynamic(1); // 0 - Explicitly disable dynamic teams, 1 - Enable dynamic teams (Faster and better optimized)
    omp_set_num_threads(numThreads); // Set the number of threads

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

    // Create the fields phi[NFields][N][N][N] and the AGD auxiliar field
    double ****phi;
    double ****psi;
    double ****phiPrev;
    phi = new double ***[NFields]; // The field is phi = sigma + pi1 + pi2 + pi3
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
    std::ifstream file0;
    std::ifstream file1;
    std::ifstream file2;
    std::ifstream file3;

    file0.open("sigma_init.dat");
    file1.open("pi1_init.dat");
    file2.open("pi2_init.dat");
    file3.open("pi3_init.dat");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                file0 >> phi[0][i][j][k];
                file1 >> phi[1][i][j][k];
                file2 >> phi[2][i][j][k];
                file3 >> phi[3][i][j][k];

                phiPrev[0][i][j][k] = phi[0][i][j][k];
                phiPrev[1][i][j][k] = phi[1][i][j][k];
                phiPrev[2][i][j][k] = phi[2][i][j][k];
                phiPrev[3][i][j][k] = phi[3][i][j][k];

                psi[0][i][j][k] = phi[0][i][j][k];
                psi[1][i][j][k] = phi[1][i][j][k];
                psi[2][i][j][k] = phi[2][i][j][k];
                psi[3][i][j][k] = phi[3][i][j][k];
            }
        }
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();
    std::cout << "Ansatz constructed" << std::endl;

    // Define the numerical quantities
    double B = 0.0;
    double EPrev = 0.0;
    double ENew = 0.0;
    double Vc = 0.0;
    double *dB = new double [3];

    // Start the minimization
    double counter = 0;
    int iter_end = 1000;
    for(int iter = 0; iter < iter_end; iter++){
        if(iter % 100 == 0){
            std::cout << "iteration = " << iter << " \t" << "velocity = " << iteration << std::endl;
            B = Baryon(phi);
            EPrev = Energy(phi);
            Vc = Virial(phi);
            std::cout << "B = " << B << "\t" << "E = " << EPrev << "\t" << "Virial = " << Vc << std::endl;
        }

        // Perform the minimization
        dEdPhi(phi, psi, phiPrev, dB);
        iteration += 1.0;
        ENew = Energy(phi);

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
double Baryon(double ****phi) {
    double B = 0.0;

    #pragma omp parallel for reduction(+:B)
    for (int i = 2; i < N - 2; i++) {
        for (int j = 2; j < N - 2; j++) {
            for (int k = 2; k < N - 2; k++) {
                double D[3][4] = {{0.0}};
                double tempB = 0.0;

                // Compute D array using private variables
                for (int a = 0; a < NFields; a++) {
                    D[0][a] = Dx(phi[a], i, j, k);
                    D[1][a] = Dy(phi[a], i, j, k);
                    D[2][a] = Dz(phi[a], i, j, k);
                }

                // Compute the Baryon term
                for (int a = 0; a < NFields; a++) {
                    for (int b = 0; b < NFields; b++) {
                        for (int c = 0; c < NFields; c++) {
                            for (int d = 0; d < NFields; d++) {
                                tempB += LeviCivita(a, b, c, d) * phi[a][i][j][k] * D[0][b] * D[1][c] * D[2][d];
                            }
                        }
                    }
                }
                // Accumulate the result using critical section
                #pragma omp critical
                B += tempB;
            }
        }
    }
    B *= -dx * dy * dz / (2.0 * M_PI * M_PI);
    return B;
}

double Energy(double ****phi) {
    double E2 = 0.0;
    double E4 = 0.0;
    double E6 = 0.0;
    double E0 = 0.0;

    #pragma omp parallel for reduction(+:E2,E4,E6,E0)
    for (int i = 2; i < N - 2; i++) {
        for (int j = 2; j < N - 2; j++) {
            for (int k = 2; k < N - 2; k++) {
                double D[3][4] = {{0.0}};
                // Compute D array using private variables
                for (int a = 0; a < NFields; a++) {
                    D[0][a] = Dx(phi[a], i, j, k);
                    D[1][a] = Dy(phi[a], i, j, k);
                    D[2][a] = Dz(phi[a], i, j, k);
                }
                double E6_local = pow(2.0 * M_PI * M_PI * BaryonDensity(phi, i, j, k), 2);
                double E0_local = 1.0 - phi[0][i][j][k];

                #pragma omp atomic
                E6 += E6_local;
                #pragma omp atomic
                E0 += E0_local;

                double E2_local = 0.0;
                double E4_local = 0.0;
                for (int a = 0; a < NFields; a++) {
                    E2_local += D[0][a] * D[0][a] + D[1][a] * D[1][a] + D[2][a] * D[2][a];

                    for (int b = 0; b < NFields; b++) {
                        for (int d1 = 0; d1 < 3; d1++) {
                            for (int d2 = 0; d2 < 3; d2++) {
                                E4_local += D[d1][a] * D[d1][a] * D[d2][b] * D[d2][b] - D[d1][a] * D[d2][b] * D[d1][b] * D[d2][a];
                            }
                        }
                    }
                }

                #pragma omp atomic
                E2 += E2_local;
                #pragma omp atomic
                E4 += E4_local;
            }
        }
    }
    double E = dx * dy * dz / (24.0 * M_PI * M_PI) * (E2 + 2.0 * E4 + c6 * E6 + c0 * E0);
    return E;
}

double Virial(double ****phi) {
    double E2 = 0.0;
    double E4 = 0.0;
    double E6 = 0.0;
    double E0 = 0.0;

    #pragma omp parallel for reduction(+:E2,E4,E6,E0)
    for (int i = 2; i < N - 2; i++) {
        for (int j = 2; j < N - 2; j++) {
            for (int k = 2; k < N - 2; k++) {
                double D[3][4] = {{0.0}};

                // Compute D array using private variables
                for (int a = 0; a < NFields; a++) {
                    D[0][a] = Dx(phi[a], i, j, k);
                    D[1][a] = Dy(phi[a], i, j, k);
                    D[2][a] = Dz(phi[a], i, j, k);
                }
                double E6_local = pow(2.0 * M_PI * M_PI * BaryonDensity(phi, i, j, k), 2);
                double E0_local = 1.0 - phi[0][i][j][k];

                #pragma omp atomic
                E6 += E6_local;
                #pragma omp atomic
                E0 += E0_local;

                double E2_local = 0.0;
                double E4_local = 0.0;

                for (int a = 0; a < NFields; a++) {
                    E2_local += D[0][a] * D[0][a] + D[1][a] * D[1][a] + D[2][a] * D[2][a];

                    for (int b = 0; b < NFields; b++) {
                        for (int d1 = 0; d1 < 3; d1++) {
                            for (int d2 = 0; d2 < 3; d2++) {
                                E4_local += D[d1][a] * D[d1][a] * D[d2][b] * D[d2][b] - D[d1][a] * D[d2][b] * D[d1][b] * D[d2][a];
                            }
                        }
                    }
                }

                #pragma omp atomic
                E2 += E2_local;
                #pragma omp atomic
                E4 += E4_local;
            }
        }
    }
    double numerator = 1.0 - 2.0 * E4 / (E2 - 3.0 * c6 * E6 + 3.0 * c0 * E0);
    return dx * dy * dz * fabs(numerator);
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

    #pragma omp parallel for shared(phi, psi, phiPrev, dB) private(D, DD, dEdphi, Vs_0, Vs_1, Vs_2, Lag6_1, Lag6_2, norma, normapsi)
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
                    dEdphi += -8*c6*Lag6_1*Lag6_2 + 2*c6*(Vs_0*dB[0] + Vs_1*dB[1] + Vs_2*dB[2]); // Sextic
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
}

// Derivative of the Baryon density
void diffB(double ****phi, int i, int j, int k, double *dB){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}; // D[coordinate][Field]
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

    double temp_dB0 = 0.0;
    double temp_dB1 = 0.0;
    double temp_dB2 = 0.0;

    #pragma omp parallel for reduction(+:temp_dB0, temp_dB1, temp_dB2)
    for(int a = 0; a<NFields; a++){
        for(int b = 0; b<NFields; b++){
            for(int c = 0; c<NFields; c++){
                for(int d = 0; d<NFields; d++){
                    double temp0 = LeviCivita(a, b, c, d) * phi[a][i][j][k] * (DD[0][0][b] * D[1][c] * D[2][d] + DD[0][1][c] * D[0][b] * D[2][d] + DD[0][2][d] * D[0][b] * D[1][c]);
                    double temp1 = LeviCivita(a, b, c, d) * phi[a][i][j][k] * (DD[1][0][b] * D[1][c] * D[2][d] + DD[1][1][c] * D[0][b] * D[2][d] + DD[1][2][d] * D[0][b] * D[1][c]);
                    double temp2 = LeviCivita(a, b, c, d) * phi[a][i][j][k] * (DD[2][0][b] * D[1][c] * D[2][d] + DD[2][1][c] * D[0][b] * D[2][d] + DD[2][2][d] * D[0][b] * D[1][c]);

                    #pragma omp atomic
                    temp_dB0 += temp0;
                    #pragma omp atomic
                    temp_dB1 += temp1;
                    #pragma omp atomic
                    temp_dB2 += temp2;
                }
            }
        }
    }
    // Fill the dB array
    dB[0] += temp_dB0;
    dB[1] += temp_dB1;
    dB[2] += temp_dB2;
}

double BaryonDensity(double ****phi, int i, int j, int k){
    double D[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}; // D[coordinate][Field]
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
