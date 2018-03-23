//
// Created by Cooper on 3/19/18.
//

#ifndef SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H
#define SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H

#include <iostream>
#include <complex>
#include <cmath>
#include "Data.h"
//#include <InvertMatrix.h>
#include <iostream>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;
namespace ublas = boost::numeric::ublas;

class Calculation {
    Data data;
private:
    //  Data data;
    int number_of_spins, lat_s, omega_length;
    complex<float> ****GF,****GF_inversive, ****GF0_inversive, ****GF0, ****GF_final_inversive, ****GF_final, ****Sigma;
    complex<float> One;
public:
    Calculation(const Data &data1) {

        data = data1;
        //construct_data_set(local_coulomb, nonlocal_exchange,
//                           chemical_potential, hopping,
//                           inversive_temperature, number_of_frequencies,
//                           sites);

        number_of_spins = data.get_number_of_spins();
        lat_s = data.get_number_of_sites();
        omega_length = data.get_number_of_freq();

        GF              = new complex<float> ***[number_of_spins];
        GF_inversive        = new complex<float> ***[number_of_spins];
        GF0              = new complex<float> ***[number_of_spins];
        GF0_inversive   = new complex<float> ***[number_of_spins];
        Sigma           = new complex<float> ***[number_of_spins];
        GF_final_inversive        = new complex<float> ***[number_of_spins];
        GF_final        = new complex<float> ***[number_of_spins];

        for (int spin = 0; spin < number_of_spins; spin++) {
            GF[spin] = new complex<float> **[lat_s];
            GF_inversive[spin] = new complex<float> **[lat_s];
            GF0[spin] = new complex<float> **[lat_s];
            GF0_inversive[spin] = new complex<float> **[lat_s];
            Sigma[spin] = new complex<float> **[lat_s];
            GF_final_inversive[spin] = new complex<float> **[lat_s];
            GF_final[spin] = new complex<float> **[lat_s];

            for (int i = 0; i < lat_s; i++) {
                GF[spin][i] = new complex<float> *[lat_s];
                GF_inversive[spin][i] = new complex<float> *[lat_s];
                GF0[spin][i] = new complex<float> *[lat_s];
                GF0_inversive[spin][i] = new complex<float> *[lat_s];
                Sigma[spin][i] = new complex<float> *[lat_s];
                GF_final_inversive[spin][i] = new complex<float> *[lat_s];
                GF_final[spin][i] = new complex<float> *[lat_s];

                for (int j = 0; j < lat_s; j++) {
                    GF[spin][i][j] = new complex<float>[omega_length];
                    GF_inversive[spin][i][j] = new complex<float>[omega_length];
                    GF0[spin][i][j] = new complex<float>[omega_length];
                    GF0_inversive[spin][i][j] = new complex<float>[omega_length];
                    Sigma[spin][i][j] = new complex<float>[omega_length];
                    GF_final_inversive[spin][i][j] = new complex<float>[omega_length];
                    GF_final[spin][i][j] = new complex<float>[omega_length];

                    for (int freq = 0; freq < omega_length; freq++) {
                        GF[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        GF_inversive[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        GF0[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        GF0_inversive[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        Sigma[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        GF_final_inversive[spin][i][j][freq] = complex<float>(0.0, 0.0);
                        GF_final[spin][i][j][freq] = complex<float>(0.0, 0.0);
                    }
                }
            }
        }
    }

    complex<float> ****get_GF() const {
        return GF;
    }

    complex<float> ****get_GF0() const {
        return GF0_inversive;
    }

    complex<float> ****get_Sigma() const {
        return Sigma;
    }

    complex<float> ****get_GF_final() const {
        return GF_final_inversive;
    }

    void construct_initial_lattice_function() {
        One = complex<float>(1.0, 0.0);
        complex<float> temp;
        GF0_inversive = new complex<float> ***[number_of_spins];
        for (int spin = 0; spin < number_of_spins; spin++) {
            GF0_inversive[spin] = new complex<float> **[lat_s];
            for (int i = 0; i < lat_s; i++) {
                GF0_inversive[spin][i] = new complex<float> *[lat_s];
                for (int j = 0; j < lat_s; j++) {
                    GF0_inversive[spin][i][j] = new complex<float>[omega_length];
                    for (int freq = 0; freq < omega_length; freq++) {
                        if (i == j) {
                            temp = data.get_omega()[freq] + data.get_chemical_potential();
                        } else {
                            temp = complex<float> (0.0, 0.0);
                        }
                        GF0_inversive[spin][i][j][freq] = temp - data.get_t_matrix()[i][j];
                    }
                }
            }
        }
    }


    void compare_GF(){
        cout << GF[0][0][0][0].imag() << " vs " << GF_final[0][0][0][0].imag() << endl;
    }

    bool test_convergency() {
        compare_GF();
        if (abs(GF[0][0][0][0].imag() - GF_final[0][0][0][0].imag()) < 0.001) {
            return true;
        } else {
            return false;
        }
    }

    void GF_takes_GO() {
        GF = GF0;
        GF_inversive = GF0_inversive;
    }

    void compute_Sigma() {
        int lat_s, len_ferm_freq, len_bos_freq, number_of_spins;
        int new_freq, new_freq_2;
        complex<float> GF_tail_kl, GF_tail_qp, GF_tail_km, GF_left_shift;
        GF_tail_kl      = complex<float> (0.0, 0.0);
        GF_tail_qp      = complex<float> (0.0, 0.0);
        GF_tail_km      = complex<float> (0.0, 0.0);
        GF_left_shift   = complex<float> (0.0, 0.0);
        lat_s = data.get_number_of_sites();
        len_ferm_freq = data.get_number_of_freq();
        len_bos_freq = data.get_number_of_freq();
        number_of_spins = data.get_number_of_spins();

        float beta;
        beta = data.get_beta();

        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int omega_n = 0; omega_n < len_ferm_freq; omega_n++) {

                        /* I term */
                        for (int spin_prime = 0; spin_prime < number_of_spins; spin_prime++) {
                            for (int k = 0; k < lat_s; k++) {
                                for (int l = 0; l < lat_s; l++) {
                                    for (int omega_n_prime = 0; omega_n_prime < len_ferm_freq; omega_n_prime++) {
                                        Sigma[spin][i][j][omega_n] +=
                                                data.get_U_matrix()[i][j][k][l] * GF[spin_prime][k][l][omega_n_prime];
                                    }
                                }
                            }
                        }

                        /* II term */

                        for (int k = 0; k < lat_s; k++) {
                            for (int l = 0; l < lat_s; l++) {
                                for (int Omega_n = 0; Omega_n < len_bos_freq; Omega_n++) {
                                    new_freq = omega_n + Omega_n;
                                    if (new_freq < len_ferm_freq) {
                                        GF_tail_kl = GF[spin][k][l][omega_n + Omega_n];
                                    } else {
                                        GF_tail_kl = complex<float>(1.0, 0.0) / data.fermionic_matsubara(new_freq);
                                    }
                                    Sigma[spin][i][j][omega_n] -= data.get_U_matrix()[i][k][j][l] * GF_tail_kl;
                                    GF_tail_kl = complex<float> (0.0, 0.0);
                                }
                            }
                        }


                        /* III term */

                        for (int spin_prime = 0; spin_prime < number_of_spins; spin_prime++) {
                            for (int n = 0; n < lat_s; n++) {
                                for (int k = 0; k < lat_s; k++) {
                                    for (int l = 0; l < lat_s; l++) {
                                        for (int m = 0; m < lat_s; m++) {
                                            for (int p = 0; p < lat_s; p++) {
                                                for (int q = 0; q < lat_s; q++) {
                                                    for (int omega_n_prime = 0;
                                                         omega_n_prime < len_ferm_freq; omega_n_prime++) {
                                                        for (int Omega_n = 0; Omega_n < len_bos_freq; Omega_n++) {
                                                            new_freq = omega_n + Omega_n;
                                                            if (new_freq < len_ferm_freq) {
                                                                GF_tail_kl = GF[spin][k][l][omega_n + Omega_n];
                                                            } else {
                                                                GF_tail_kl = complex<float>(1.0, 0.0) /
                                                                             data.fermionic_matsubara(new_freq);
                                                            }
                                                            new_freq = omega_n_prime + Omega_n;
                                                            if (new_freq < len_ferm_freq) {
                                                                GF_tail_qp = GF[spin_prime][q][p][new_freq];
                                                            } else {
                                                                GF_tail_qp = complex<float>(1.0, 0.0) /
                                                                             data.fermionic_matsubara(new_freq);
                                                            }
                                                            Sigma[spin][i][j][omega_n] +=
                                                                    data.get_U_matrix()[i][k][n][p] *
                                                                    data.get_U_matrix()[l][j][q][m] *
                                                                    GF_tail_kl * GF_tail_qp *
                                                                    GF[spin_prime][n][m][omega_n_prime] /
                                                                    complex<float>(pow(beta, 2), 0.0);

                                                            GF_tail_kl = complex<float> (0.0, 0.0);
                                                            GF_tail_qp = complex<float> (0.0, 0.0);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }


                        /* IV term */

                        for (int n = 0; n < lat_s; n++) {
                            for (int k = 0; k < lat_s; k++) {
                                for (int l = 0; l < lat_s; l++) {
                                    for (int m = 0; m < lat_s; m++) {
                                        for (int p = 0; p < lat_s; p++) {
                                            for (int q = 0; q < lat_s; q++) {
                                                for (int Omega_n = 0; Omega_n < len_bos_freq; Omega_n++) {
                                                    new_freq = omega_n + Omega_n;
                                                    if (new_freq < len_ferm_freq) {
                                                        GF_tail_km = GF[spin][k][m][omega_n + Omega_n];
                                                    } else {
                                                        GF_tail_km = complex<float>(1.0, 0.0) /
                                                                     data.fermionic_matsubara(new_freq);
                                                    }
                                                    new_freq_2 = omega_n - Omega_n;
                                                    if (new_freq_2 > 0) {
                                                        GF_left_shift = GF[spin][n][l][omega_n - Omega_n];
                                                    } else {
                                                        GF_left_shift = GF[spin][n][l][-(omega_n - Omega_n)];
                                                    }
                                                    Sigma[spin][i][j][omega_n] -= data.get_U_matrix()[i][k][n][p] *
                                                                                  data.get_U_matrix()[l][j][q][m] *
                                                                                  GF_tail_km * GF_left_shift *
                                                                                  GF[spin][q][p][omega_n] /
                                                                                  complex<float>(beta, 0.0);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

    }

    void Dyson_equation() {
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        GF_final_inversive[spin][i][j][freq] =  GF_inversive[spin][i][j][freq] - Sigma[spin][i][j][freq];
                    }
                }
            }
        }
    }

    void GF_takes_GF_final() {
        GF              = GF_final;
        GF_inversive    = GF_final_inversive;
    }

    void print_Sigma_in_file() {
        ofstream out;
        for (int i = 0; i < lat_s; i++) {
            for (int j = 0; j < lat_s; j++) {
                out.open("Sigma_" + to_string(i) + "_" + to_string(j) + ".txt");
                if (out.is_open()) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        out << data.get_omega()[freq].imag() << " " << Sigma[0][i][j][freq].real() << " "
                            << Sigma[0][i][j][freq].imag() << " ";
                        out << Sigma[1][i][j][freq].real() << " " << Sigma[1][i][j][freq].imag() << endl;
                    }
                }
                out.close();
            }
        }
    }

    void print_GF_in_file() {
        ofstream out;
        for (int i = 0; i < lat_s; i++) {
            for (int j = 0; j < lat_s; j++) {
                out.open("GF_" + to_string(i) + "_" + to_string(j) + ".txt");
                if (out.is_open()) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        out << data.get_omega()[freq].imag() << " " << GF[0][i][j][freq].real() << " "
                            << GF[0][i][j][freq].imag() << " ";
                        out << GF[1][i][j][freq].real() << " " << GF[1][i][j][freq].imag() << endl;
                    }
                }
                out.close();
            }
        }
    }

    void print_GF0_in_file(){
        ofstream out;

            for (int i = 0; i < lat_s; i++){
                for(int j = 0; j < lat_s; j++){
                    out.open("GF0_" + to_string(i) + "_" + to_string(j) + ".txt");
                    if (out.is_open()) {
                        for (int freq = 0; freq < omega_length; freq++) {
                            out << data.get_omega()[freq].imag() << " " << GF0[0][i][j][freq].real() << " "
                                << GF0[0][i][j][freq].imag() << " ";
                            out << GF0[1][i][j][freq].real() << " " << GF0[1][i][j][freq].imag() << endl;
                        }
                    }
                    out.close();
                }
            }

    }

    void clean_memory_GF(){

        for(int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++){
                for(int j = 0; j < lat_s; j++){
                    delete[] GF[spin][i][j];
                    delete[] GF_inversive[spin][i][j];
                    delete[] GF0[spin][i][j];
                    delete[] GF0_inversive[spin][i][j];
                    delete[] GF_final_inversive[spin][i][j];
                    delete[] GF_final[spin][i][j];
                    delete[] Sigma[spin][i][j];
                }
                delete[] GF[spin][i];
                delete[] GF_inversive[spin][i];
                delete[] GF0[spin][i];
                delete[] GF0_inversive[spin][i];
                delete[] GF_final_inversive[spin][i];
                delete[] GF_final[spin][i];
                delete[] Sigma[spin][i];
            }
            delete[] GF[spin];
            delete[] GF_inversive[spin];
            delete[] GF0[spin];
            delete[] GF0_inversive[spin];
            delete[] GF_final_inversive[spin];
            delete[] GF_final[spin];
            delete[] Sigma[spin];
        }
        delete[] GF;
        delete[] GF_inversive;
        delete[] GF0;
        delete[] GF0_inversive;
        delete[] GF_final_inversive;
        delete[] GF_final;
        delete[] Sigma;
    }

    template<class T>
    bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        // create a working copy of the input
        matrix<T> A(input);
        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());
        // perform LU-factorization
        int res = lu_factorize(A,pm);
        if( res != 0 )
            return false;
        // create identity matrix of "inverse"
        inverse.assign(ublas::identity_matrix<T>(A.size1()));
        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);
        return true;
    }

    void inverse_matrix(){
        int size = 2;
        bool inverted;

        matrix<float> start_matrix(size, size);

        start_matrix(0, 0) = 3.0;
        start_matrix(0, 1) = 4.0;
        start_matrix(1, 0) = 5.0;
        start_matrix(1, 1) = 6.0;

        matrix<float> inversion(size, size);

        inverted = InvertMatrix(start_matrix, inversion);
        cout << start_matrix << endl;
        cout << inversion << endl;

    }

    void inverse_matrix_to_G0(){
        int size = 2;
        bool inverted;

        matrix<complex<float>> start_matrix(lat_s, lat_s);
        matrix<complex<float>> inversion(lat_s, lat_s);

        for(int spin = 0; spin < number_of_spins; spin++){
            for(int freq = 0; freq < omega_length; freq++){
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        start_matrix(i, j) = GF0_inversive[spin][i][j][freq];
                    }
                }
                inverted = InvertMatrix(start_matrix, inversion);
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        GF0[spin][i][j][freq] = inversion(i, j);
                    }
                }
            }
        }

        inverted = InvertMatrix(start_matrix, inversion);
//        cout << start_matrix << endl;
//        cout << inversion << endl;

    }

    void inverse_matrix_to_GF_final(){
        int size = 2;
        bool inverted;

        matrix<complex<float>> start_matrix(lat_s, lat_s);
        matrix<complex<float>> inversion(lat_s, lat_s);

        for(int spin = 0; spin < number_of_spins; spin++){
            for(int freq = 0; freq < omega_length; freq++){
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        start_matrix(i, j) = GF_final_inversive[spin][i][j][freq];
                    }
                }
                inverted = InvertMatrix(start_matrix, inversion);
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        GF_final[spin][i][j][freq] = inversion(i, j);
                    }
                }
            }
        }

        inverted = InvertMatrix(start_matrix, inversion);
//        cout << start_matrix << endl;
//        cout << inversion << endl;

    }

    void inverse_matrix_GF(){
        int size = 2;
        bool inverted;

        matrix<complex<float>> start_matrix(lat_s, lat_s);
        matrix<complex<float>> inversion(lat_s, lat_s);

        for(int spin = 0; spin < number_of_spins; spin++){
            for(int freq = 0; freq < omega_length; freq++){
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        start_matrix(i, j) = GF[spin][i][j][freq];
                    }
                }
                inverted = InvertMatrix(start_matrix, inversion);
                for(int i = 0; i < lat_s; i++){
                    for(int j = 0; j < lat_s; j++){
                        GF_inversive[spin][i][j][freq] = inversion(i, j);
                    }
                }
            }
        }

        inverted = InvertMatrix(start_matrix, inversion);
//        cout << start_matrix << endl;
//        cout << inversion << endl;

    }
};


#endif //SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H
