#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <fstream> 
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace chrono;

double upperBound = 8.0;
int N_map[] = {65,129};
int N = 65;
const double dx = upperBound / ((double)N - 1);
double mu_inverse = (2.0) / (1836.152673);
const double D = 0.1744;
const double R_eq = 1.40201;
const double beta = 1.02764;
double dx_2_inverse = (1.0 / (dx*dx));

double V(double r);
double renormalized_hamilton_matrix(int i, int j);

//....oO00o........oO00o........oO00o........oO00o........oO00o........oO00o.... 

int main() {

    ofstream eigenvalues_out("eigenvalues.txt");
    for(int k = 0; k < 2; k++) {
        N = N_map[k]; 
        /** ********************** a) "simple" derivative ***********************
         *  b ... diagonal of the tridiagonal hamilton matrix
         *  a ... sub diagonal of the tridiagonal hamilton matrix
        */

        auto start = chrono::high_resolution_clock::now();

        VectorXd b(N);
        VectorXd a(N-1);

        b[0] = mu_inverse * dx_2_inverse + V(0.0);

        double* V_test = new double[N];
        for(int i = 1; i < N; i++){
            b[i] = mu_inverse * dx_2_inverse + V((i)*dx);
            a[i-1] = -0.5 * mu_inverse * dx_2_inverse;
        }
        
        SelfAdjointEigenSolver<MatrixXd> hamiltonSolver_a;

        hamiltonSolver_a.computeFromTridiagonal(b,a);
        
        auto stop = high_resolution_clock::now();
        auto duration_a = duration_cast<milliseconds>(stop - start);

        /** ********************** b) Fourier grid method **********************
        */

        start = chrono::high_resolution_clock::now();

        MatrixXd H_0(N, N);

        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                H_0(i,j) = renormalized_hamilton_matrix(i,j);
            }
        }

        SelfAdjointEigenSolver<MatrixXd> hamiltonSolver_b;

        hamiltonSolver_b.compute(H_0);

        stop = chrono::high_resolution_clock::now();
        auto duration_b = duration_cast<milliseconds>(stop - start);

        /** ********************** Output ***********************
         * Output file: eigenvalues.txt
         * 1. column: first 10 eigenvalues from a)
         * 2. column: relative difference to values from Table II (J. Chem. Phys. 91, 3571 (1989))
         * 3. column: first 10 eigenvalues from b)
         * 4. column: relative difference to values from Table II (J. Chem. Phys. 91, 3571 (1989))
         * 5. column: relative difference between mean value of a) and b) to the analytical values from Table II (J. Chem. Phys. 91, 3571 (1989))
         * 6. column: relative difference between values a and b
        */

        string filename_computed = "n" + to_string(N) + "_eigenvalues.txt";
        string filename_analytic = "analytical_eigenvalues.txt";
        ifstream inputFile_com(filename_computed);
        ifstream inputFile_ana(filename_analytic);
        string line;
        double in_value;
        
        if (!inputFile_com.is_open()) {
            cout << "Failed to open " << filename_computed << " file" << endl;
            break;
        } else if(!inputFile_ana.is_open()) {
            cout << "Failed to open " << filename_analytic << " file" << endl;
            break;
        }
        double diff;

        eigenvalues_out << "N = " << N << ":" << endl;
        for (int i = 0; i < 10; i++){
            double eigenvalue_a = hamiltonSolver_a.eigenvalues()(i);
            double eigenvalue_b = hamiltonSolver_b.eigenvalues()(i);
            if(getline(inputFile_com, line)) {
                in_value = stod(line);
            }
            
            eigenvalues_out << hamiltonSolver_a.eigenvalues()(i) << "\t\t";
            diff = (eigenvalue_a - in_value)/in_value;
            eigenvalues_out << ((diff > 0.0) ? "+ " : "") << diff << "\t\t"
                            << hamiltonSolver_b.eigenvalues()(i) << "\t\t";
            diff = (eigenvalue_b - in_value)/in_value;
            eigenvalues_out << ((diff > 0.0) ? "+ " : "") << diff << "\t\t";
            if(getline(inputFile_ana, line)) {
                in_value = stod(line);
            }
            diff = (((eigenvalue_a + eigenvalue_b) / 2.0) - in_value)/in_value;
            eigenvalues_out << ((diff > 0.0) ? "+ " : "") << diff << "\t\t"
                            << (eigenvalue_b - eigenvalue_a)/eigenvalue_a << endl;
        }
        eigenvalues_out << "Performance a: " << duration_a.count() << " ms" << endl 
                        << "Performance b [ms]: " << duration_b.count() << " ms" << endl;
        eigenvalues_out << endl << endl;

        inputFile_com.close();

    }
    eigenvalues_out.close();
    
    return 0;

}

//....oO00o........oO00o........oO00o........oO00o........oO00o........oO00o.... 

double V(double r) {
    return D * pow((1.0 - exp(-beta*(r - R_eq))),2); 
}

//....oO00o........oO00o........oO00o........oO00o........oO00o........oO00o.... 

double renormalized_hamilton_matrix(int i, int j) {
    double H_0=0.0;
    int n = (int)((N-1)/2);
    double N_inverse = 1.0/N;
    double dx_inverese = 1.0 / dx;

    for(int l=0; l < n; l++){
        H_0 += cos(l*2*M_PI*(i-j)*N_inverse)*2*mu_inverse*(M_PI*M_PI*l*l*N_inverse*N_inverse*dx_inverese*dx_inverese);
    } 
    
    H_0 = H_0 * 2.0 * N_inverse;
    if(i == j)
        H_0 += V((double)i*dx);

    return H_0;
}
