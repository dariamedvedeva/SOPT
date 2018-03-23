//
// Created by Cooper on 3/19/18.
//
#include <map>
#include <complex>
#include <vector>
#include <math.h>
#include <cmath>

#ifndef SECON_ORGER_DIAGRAM_SIGMA_DATA_H
#define SECON_ORGER_DIAGRAM_SIGMA_DATA_H

using namespace std;

class Data {
private:
    float U, J, t, mu, beta;
    int N;
    int number_of_sites, number_of_spins = 2;
    map<int, vector<int> > connections, connections_10_10;
    complex<float> **t_matrix, *omega, *nu, ****U_matrix;

public:
    Data(){
        set_U(0.0);
        set_J(0.0);
        set_chemical_potential(0.0);
        set_hopping(0.0);
        set_beta(0.0);
        set_number_of_freq(0);
        set_number_of_sites(0);
    }

    void init_parameters(float local_coulomb, float nonlocal_exchange, float chemical_potential, float hopping,
                         float inversive_temperature, int number_of_frequencies, int sites){
        set_U(local_coulomb);
        set_J(nonlocal_exchange);
        set_chemical_potential(chemical_potential);
        set_hopping(hopping);
        set_beta(inversive_temperature);
        set_number_of_freq(number_of_frequencies);
        set_number_of_sites(sites);

    }

    void print_data(){
        cout << "/n Parameters: "                   << get_U() << endl;
        cout << "Local Coulomb parameter = "        << get_U() << endl;
        cout << "Non-local exchange parameter = "   << get_J() << endl;
        cout << "Chemical potential = "             << get_chemical_potential() << endl;
        cout << "Beta = "                           << get_beta() << endl;
        cout << "The number of frequencies = "      << get_number_of_freq() << endl;
        cout << "The number of sites = "            << get_number_of_sites() << endl;
        cout << endl;
    }

    void construct_hopping_matrix(){
        t_matrix = new complex <float> *[get_number_of_sites()];
        for(int i = 0; i < get_number_of_sites(); i++) {
            t_matrix[i] = new complex<float>[get_number_of_sites()];
        }
        for (int i = 0; i < get_number_of_sites(); i++){
            for(int j = 0; j < get_number_of_sites(); j++){
                if (i != j){
                    for (int p = 0; p < sqrt(get_number_of_sites()); p++) {
                        if (j == connections.at(i).at(p)) {
                            t_matrix[i][j] = complex<float> (get_t(), 0.0);
                        }
                    }
                }
            }
        }
    }

    void print_hopping_matrix(){
        cout << "Hopping matrix:" << endl;
        for (int i = 0; i < get_number_of_sites(); i++){
            for(int j = 0; j < get_number_of_sites(); j++){
                cout << t_matrix[i][j] << '\t';
            }
            cout << endl;
        }
        cout << endl;
    }

    void construct_connections(){
        /* It's a numeric magic, dude) */
        vector<int> con1 = {1, 2};
        vector<int> con2 = {0, 3};
        vector<int> con3 = {0, 3};
        vector<int> con4 = {1, 2};
        connections =  { {0, con1}, {1, con2}, {2, con3}, {3, con4} };
    }

    void construct_connection_10_10(){
        /* Lattice 10 x 10*/
        
    }

    void clear_memory_t_matrix(){
        for(int i; i < get_number_of_sites(); i++) {
            delete[] t_matrix[i];
        }
        delete[] t_matrix;
    }

    void frequencies(){
        omega   = new complex<float> [get_number_of_freq()];
        nu      = new complex<float> [get_number_of_freq()];
        for(int i = 0; i < get_number_of_freq(); i++){
            omega[i]    = fermionic_matsubara(i);
            nu[i]       = bosonic_matsubara(i);
        }
    }

    void construct_U_matrix() {
        int lat_s = get_number_of_sites();

        U_matrix = new complex<float> ***[lat_s];
        for (int i = 0; i < lat_s; i++) {
            U_matrix[i] = new complex<float> **[lat_s];
            for (int j = 0; j < lat_s; j++) {
                U_matrix[i][j] = new complex<float> *[lat_s];
                for (int k = 0; k < lat_s; k++) {
                    U_matrix[i][j][k] = new complex<float>[lat_s];
                    for (int l = 0; l < lat_s; l++) {
                        U_matrix[i][j][k][l] = complex<float>(get_J(), 0.0);
                        if ((i == k) & (j == l) & (i != j)) {
                            for (int p = 0; p < sqrt(get_number_of_sites()); p++) {
                                if (j == connections.at(i).at(p)) {
                                    U_matrix[i][j][k][l] = complex<float>(get_J(), 0.0);
                                    /* Works */
                                    //cout << "Nonlocal: [" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                                }
                            }
                        }
                        if ((i == j) & (i == k) & (i == l)) {
                            U_matrix[i][j][k][l] = complex<float>(U / 2.0, 0.0);
                            /* Works */
                            //cout << "Local: [" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                        }
                    }
                }
            }
        }
    }

    void print_U_matrix(){
        int lat_s = get_number_of_sites();
        cout << "U matrix:" << endl;
        for(int i = 0; i < lat_s; i++){
            for (int j = 0; j < lat_s; j++){
                for(int k = 0; k < lat_s; k++){
                    for(int l = 0; l < lat_s; l++){
                        if ((i == j) & (i == k) & (i == l)){
                            cout << "Local: \t\t[" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                        }
                        if ((i == k) & (j == l) & (i != j)){
                            for (int p = 0; p < sqrt(get_number_of_sites()); p++) {
                                if (j == connections.at(i).at(p)) {
                                    cout << "Nonlocal: \t[" << i << "][" << j << "][" << k << "][" << l << "] = "
                                         << U_matrix[i][j][k][l] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        cout << endl;
    }

    complex<float> fermionic_matsubara(int n){
        return float(M_PI * (2 * n + 1) / get_beta()) * complex<float> (0, 1);
    }

    complex<float> bosonic_matsubara(int n){
        return float(M_PI * 2 * n / get_beta()) * complex<float> (0, 1);
    }

    float get_U() const {
        return U;
    }

    float get_J() const {
        return J;
    }

    float get_t() const {
        return t;
    }

    complex<float> get_chemical_potential() const {
        return complex<float> (mu, 0.0);
    }

    float get_beta() const {
        return beta;
    }

    int get_number_of_freq() const {
        return N;
    }

    int get_number_of_sites() const {
        return number_of_sites;
    }

    int get_number_of_spins() const {
        return number_of_spins;
    }

    complex<float> ****get_U_matrix() const {
        return U_matrix;
    }

    complex<float> **get_t_matrix() const {
        return t_matrix;
    }

    complex<float> *get_omega() const {
        return omega;
    }

    complex<float> *get_nu() const {
        return nu;
    }

private:

    void set_U(float local_coulomb){
        U = local_coulomb;
    }

    void set_J(float nonlocal_exchange){
        J = nonlocal_exchange;
    }

    void set_chemical_potential(float chemical_potential){
        mu = chemical_potential;
    }

    void set_hopping(float hopping){
        t = hopping;
    }

    void set_beta(float inversive_temperature){
        beta = inversive_temperature;
    }

    void set_number_of_freq(int number_of_frequencies){
        N = number_of_frequencies;
    }

    void set_number_of_sites(int sites){
        number_of_sites = sites;
    }
};
#endif //SECON_ORGER_DIAGRAM_SIGMA_DATA_H