//
// Created by Cooper on 3/19/18.
//

#ifndef SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H
#define SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H

#include <iostream>
#include <complex>
#include <math.h>
#include <cmath>
#include "Data.h"
//#include <InvertMatrix.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "boost/multi_array.hpp"
#include <cassert>

using namespace std;
//using namespace boost::numeric::ublas;
namespace ublas = boost::numeric::ublas;

class Calculation {
    //Data data;
private:
    Data data;
    int number_of_spins, lat_s, omega_length;

    complex<double> One;

    vector< vector < vector < vector <complex<double> > > > > GF;
    vector< vector < vector < vector <complex<double> > > > > GF0;
    vector< vector < vector < vector <complex<double> > > > > GF_final;
    vector< vector < vector < vector <complex<double> > > > > Sigma;
    vector< vector < vector < vector <complex<double> > > > > GF_inversive;
    vector< vector < vector < vector <complex<double> > > > > GF0_inversive;
    vector< vector < vector < vector <complex<double> > > > > GF_final_inversive;

public:
    Calculation(const Data &data1) {

        data = data1;

        number_of_spins = data.get_number_of_spins();
        lat_s = data.get_number_of_sites();
        omega_length = data.get_number_of_freq();

        cout << number_of_spins << ", " << lat_s << ", " << omega_length << endl;

        resize_functions(number_of_spins, lat_s, omega_length);
        cout << GF0.size() << endl;

        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        GF[spin][i][j][freq]            = std::complex <double > (0.0, 0.0);
                        GF_inversive[spin][i][j][freq]  = std::complex <double >(0.0, 0.0);
                        GF0[spin][i][j][freq]           = std::complex <double >(0.0, 0.0);
                        GF0_inversive[spin][i][j][freq] = std::complex <double >(0.0, 0.0);
                        Sigma[spin][i][j][freq]         = std::complex <double >(0.0, 0.0);
                        GF_final_inversive[spin][i][j][freq] = std::complex <double >(0.0, 0.0);
                        GF_final[spin][i][j][freq]      = std::complex <double >(0.0, 0.0);
                    }
                }
            }
        }
    }

    void resize_functions(int spins, int lattice_size, int frequencies){

        GF0.resize(number_of_spins);
        for(auto &lvl_1 : GF0){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        GF.resize(number_of_spins);
        for(auto &lvl_1 : GF){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        GF_final.resize(number_of_spins);
        for(auto &lvl_1 : GF_final){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        Sigma.resize(number_of_spins);
        for(auto &lvl_1 : Sigma){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        GF0_inversive.resize(number_of_spins);
        for(auto &lvl_1 : GF0_inversive){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        GF_inversive.resize(number_of_spins);
        for(auto &lvl_1 : GF_inversive){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }

        GF_final_inversive.resize(number_of_spins);
        for(auto &lvl_1 : GF_final_inversive){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(omega_length);
                }
            }
        }
    }

    vector < vector < vector < vector <complex<double> > > > > get_GF() const {
        return GF;
    }

    vector < vector < vector < vector <complex<double> > > > > get_GF0_inversive() const {
        return GF0_inversive;
    }

    vector < vector < vector < vector <complex<double> > > > > get_Sigma() const {
        return Sigma;
    }

    vector < vector < vector < vector <complex<double> > > > > get_GF_final_inversive() const {
        return GF_final_inversive;
    }

    void construct_initial_lattice_function() {
        One = complex<double>(1.0, 0.0);
        complex<double> temp;
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        if (i == j) {
                            temp = data.get_omega()[freq] + data.get_chemical_potential();
                        } else {
                            temp = complex<double> (0.0, 0.0);
                        }
                        GF0_inversive[spin][i][j][freq] = temp - data.get_t_matrix()[i][j];
                    }
                }
            }
        }
        inverse_matrix_to_G0();
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
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        GF[spin][i][j][freq]            = GF0[spin][i][j][freq];
                        GF_inversive[spin][i][j][freq]  = GF0_inversive[spin][i][j][freq];
                    }
                }
            }
        }
        set_to_zero_GF_final();
    }

    void compute_Sigma() {
        int lat_s, len_ferm_freq, len_bos_freq, number_of_spins;
        int new_freq, new_freq_2;
        complex<double> GF_tail_kl, GF_tail_qp, GF_tail_km, GF_left_shift;
        GF_tail_kl      = complex<double>(0.0, 0.0);
        GF_tail_qp      = complex<double>(0.0, 0.0);
        GF_tail_km      = complex<double>(0.0, 0.0);
        GF_left_shift    = complex<double>(0.0, 0.0);
        lat_s           = data.get_number_of_sites();
        len_ferm_freq   = data.get_number_of_freq();
        len_bos_freq    = data.get_number_of_freq();
        number_of_spins = data.get_number_of_spins();

        vector < vector < vector < vector < std::complex < double > > > > > Um;

        Um.resize(lat_s);
        for(auto &lvl_1 : Um){
            lvl_1.resize(lat_s);
            for(auto &lvl_2 : lvl_1){
                lvl_2.resize(lat_s);
                for(auto &lvl_3 : lvl_2){
                    lvl_3.resize(lat_s);
                }
            }
        }

        for (int i = 0; i < lat_s; i++) {
            for (int j = 0; j < lat_s; j++) {
                for (int k = 0; k < lat_s; k++) {
                    for (int l = 0; l < lat_s; l++) {
                        Um[i][j][k][l] = data.get_U_matrix()[i][j][k][l];
                    }
                }
            }
        }



        double beta;
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
                                                Um[i][j][k][l] * GF[spin_prime][k][l][omega_n_prime];
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
                                        GF_tail_kl = complex<double>(1.0, 0.0) / data.fermionic_matsubara(new_freq);
                                    }
                                    Sigma[spin][i][j][omega_n] -= Um[i][k][j][l] * GF_tail_kl;
                                    GF_tail_kl = complex<double>(0.0, 0.0);
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
                                                                GF_tail_kl = complex<double>(1.0, 0.0) /
                                                                             data.fermionic_matsubara(new_freq);
                                                            }
                                                            new_freq = omega_n_prime + Omega_n;
                                                            if (new_freq < len_ferm_freq) {
                                                                GF_tail_qp = GF[spin_prime][q][p][new_freq];
                                                            } else {
                                                                GF_tail_qp = complex<double>(1.0, 0.0) /
                                                                             data.fermionic_matsubara(new_freq);
                                                            }
                                                            Sigma[spin][i][j][omega_n] +=
                                                                    Um[i][k][n][p] *
                                                                    Um[l][j][q][m] *
                                                                    GF_tail_kl * GF_tail_qp *
                                                                    GF[spin_prime][n][m][omega_n_prime] /
                                                                    complex<double>(pow(beta, 2), 0.0);

                                                            GF_tail_kl = complex<double>(0.0, 0.0);
                                                            GF_tail_qp = complex<double>(0.0, 0.0);
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
                                                        GF_tail_km = complex<double>(1.0, 0.0) /
                                                                     data.fermionic_matsubara(new_freq);
                                                    }
                                                    new_freq_2 = omega_n - Omega_n;
                                                    if (new_freq_2 > 0) {
                                                        GF_left_shift = GF[spin][n][l][omega_n - Omega_n];
                                                    } else {
                                                        GF_left_shift = GF[spin][n][l][-(omega_n - Omega_n)];
                                                    }
                                                    Sigma[spin][i][j][omega_n] -= Um[i][k][n][p] *
                                                                                  Um[l][j][q][m] *
                                                                                  GF_tail_km * GF_left_shift *
                                                                                  GF[spin][q][p][omega_n] /
                                                                                  complex<double>(beta, 0.0);
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
            cout << "Sigma spin " << spin << " is done." << endl;
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
        inverse_matrix_to_GF_final();
    }

    void GF_takes_GF_final() {
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        GF[spin][i][j][freq]            = GF_final[spin][i][j][freq];
                        GF_inversive[spin][i][j][freq]  = GF_final_inversive[spin][i][j][freq];
                    }
                }
            }
        }
        set_to_zero_GF_final();
    }

    void set_to_zero_GF_final(){
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        GF_final_inversive[spin][i][j][freq]    = complex<double>(0.0, 0.0);
                        GF_final[spin][i][j][freq]              = complex<double>(0.0, 0.0);
                    }
                }
            }
        }
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

//    void clean_memory_GF(){
//
//        for(int spin = 0; spin < number_of_spins; spin++) {
//            for (int i = 0; i < lat_s; i++){
//                for(int j = 0; j < lat_s; j++){
//                    delete[] GF[spin][i][j];
//                    delete[] GF_inversive[spin][i][j];
//                    delete[] GF0[spin][i][j];
//                    delete[] GF0_inversive[spin][i][j];
//                    delete[] GF_final_inversive[spin][i][j];
//                    delete[] GF_final[spin][i][j];
//                    delete[] Sigma[spin][i][j];
//                }
//                delete[] GF[spin][i];
//                delete[] GF_inversive[spin][i];
//                delete[] GF0[spin][i];
//                delete[] GF0_inversive[spin][i];
//                delete[] GF_final_inversive[spin][i];
//                delete[] GF_final[spin][i];
//                delete[] Sigma[spin][i];
//            }
//            delete[] GF[spin];
//            delete[] GF_inversive[spin];
//            delete[] GF0[spin];
//            delete[] GF0_inversive[spin];
//            delete[] GF_final_inversive[spin];
//            delete[] GF_final[spin];
//            delete[] Sigma[spin];
//        }
//        delete[] GF;
//        delete[] GF_inversive;
//        delete[] GF0;
//        delete[] GF0_inversive;
//        delete[] GF_final_inversive;
//        delete[] GF_final;
//        delete[] Sigma;
//    }

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

        ublas::matrix< std::complex<double> > start_matrix(size, size);
        //cout << complex<double > (3.0, 1.0);

        start_matrix(0, 0) = std::complex< double > (3.0, 1.0);
        start_matrix(0, 1) = std::complex< double > (4.0, 1.0);
        start_matrix(1, 0) = std::complex< double > (5.0, 1.0);
        start_matrix(1, 1) = std::complex< double > (6.0, 1.0);
//
//        ublas::matrix<std::complex <double>> inversion(size, size);
//
//        inverted = InvertMatrix(start_matrix, inversion);
//        cout << start_matrix << endl;
//        cout << inversion << endl;

    }

    void inverse_matrix_to_G0(){
        bool inverted;

        ublas::matrix< std::complex<double> > start_matrix(lat_s, lat_s);
        ublas::matrix< std::complex<double> > inversion(lat_s, lat_s);

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
    }

    void inverse_matrix_to_GF_final(){
        bool inverted;

        ublas::matrix< std::complex<double> > start_matrix(lat_s, lat_s);
        ublas::matrix< std::complex<double> > inversion(lat_s, lat_s);

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
    }

//    void inverse_matrix_GF(){
//        bool inverted;
//
//        matrix<complex<double>> start_matrix(lat_s, lat_s);
//        matrix<complex<double>> inversion(lat_s, lat_s);
//
//        for(int spin = 0; spin < number_of_spins; spin++){
//            for(int freq = 0; freq < omega_length; freq++){
//                for(int i = 0; i < lat_s; i++){
//                    for(int j = 0; j < lat_s; j++){
//                        start_matrix(i, j) = GF[spin][i][j][freq];
//                    }
//                }
//                inverted = InvertMatrix(start_matrix, inversion);
//                for(int i = 0; i < lat_s; i++){
//                    for(int j = 0; j < lat_s; j++){
//                        GF_inversive[spin][i][j][freq] = inversion(i, j);
//                    }
//                }
//            }
//        }
//    }
};


#endif //SECON_ORGER_DIAGRAM_SIGMA_GREEN_FUNCTION_H
