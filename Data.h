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
    double U, J, t, mu, beta;
    int N;
    int number_of_sites, number_of_spins = 2;
    map<int, vector<int> > connections;
    vector < std :: complex < double > > omega, nu;
    vector < vector < std::complex < double > > > t_matrix;
    vector < vector < vector < vector < std::complex < double > > > > > U_matrix;

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

    void init_parameters(double local_coulomb, double nonlocal_exchange, double chemical_potential, double hopping,
                         double inversive_temperature, int number_of_frequencies, int sites){
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
        int size = connections.size();
        t_matrix.resize(size);
        for(auto &lvl_1 : t_matrix){
            lvl_1.resize(size);
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                t_matrix[i][j] = std::complex<double> (0.0, 0.0);
            }
        }


        for (int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                if (i != j){
                    for (int p = 0; p < connections.at(i).size(); p++) {
                        if (j == connections.at(i).at(p)) {
                            t_matrix[i][j] = std::complex<double> (get_t(), 0.0);
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
        if (get_number_of_sites()  == 4){
            vector<int> con1 = {1, 2};
            vector<int> con2 = {0, 3};
            vector<int> con3 = {0, 3};
            vector<int> con4 = {1, 2};

            connections =  { {0, con1}, {1, con2}, {2, con3}, {3, con4} };
        } else if(get_number_of_sites()  == 9) {
            vector<int> con0 = {1, 2, 3, 6};
            vector<int> con1 = {0, 2, 4, 7};
            vector<int> con2 = {0, 1, 5, 8};
            vector<int> con3 = {0, 4, 5, 6};
            vector<int> con4 = {1, 3, 5, 7};
            vector<int> con5 = {2, 3, 4, 8};
            vector<int> con6 = {0, 3, 7, 8};
            vector<int> con7 = {1, 4, 6, 8};
            vector<int> con8 = {2, 5, 6, 7};
            connections =  { {0, con0}, {1, con1}, {2, con2}, {3, con3}, {4, con4}, {5, con5}, {6, con6}, {7, con7}, {8, con8} };
        } else if(get_number_of_sites() == 100) {
            /* It's a numeric magic, dude) */
            typedef pair <const int, vector<int>> Int_Pair;
            for(int i = 0; i < get_number_of_sites(); i++){
                vector<int> temp;
               // cout << "site: " << i << ":\t";
                if (temp.empty()) {
                    for (int j = 0; j < get_number_of_sites(); j++) {
                        /* bottom and top neighbours */
                        if ((j % 10 == i % 10) & (abs(i / 10 - j / 10) == 1)) {
                            // cout << j << " ";
                            temp.push_back(j);
                            /* neighbours from the left and from the ride sides */
                        } else if (abs(i - j) == 1 & i / 10 == j / 10) {
                            // cout << j << " ";
                            temp.push_back(j);
                            /* bottom and top lines */
                        } else if (abs(i - j) == 90) {
                            // cout << j << " ";
                            temp.push_back(j);
                            /* left and right boundaries */
                        } else if ((abs(i - j) == 9) & (i / 10 == j / 10) & (i % 10 == 0 || i % 10 == 9)) {
                            // cout << j << " ";
                            temp.push_back(j);
                        }
                    }
                    connections.insert(Int_Pair(i, temp));
                    temp.clear();
                }
            }
        }
    }

    void print_connections(){
        cout << "Print connections" << endl;
        cout << "cite" << "\t neighbours" << endl;
        for (auto it = connections.begin(); it != connections.end(); ++it) {
            cout << it->first << " : {";
            for (vector<int>::const_iterator arr = it->second.begin(); arr != it->second.end(); ++arr) {
                cout << *arr << " ";
            }
            cout << "}" << endl;
        }
    }

    void frequencies(){
        int frequencies = get_number_of_freq();
        omega.resize(frequencies, std::complex<double> (0.0, 0.0));
        nu.resize(frequencies, std::complex<double> (0.0, 0.0));

        for(int i = 0; i < get_number_of_freq(); i++){
            omega[i]    = fermionic_matsubara(i);
            nu[i]       = bosonic_matsubara(i);
        }

    }

    void construct_U_matrix() {
        int lat_s = connections.size();

        U_matrix.resize(lat_s);
        for(auto &lvl_1 : U_matrix){
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
                        U_matrix[i][j][k][l] = std::complex <double> (0.0, 0.0);
                    }
                }
            }
        }

        for (int i = 0; i < lat_s; i++) {
            for (int j = 0; j < lat_s; j++) {
                for (int k = 0; k < lat_s; k++) {
                    for (int l = 0; l < lat_s; l++) {
                        U_matrix[i][j][k][l] = complex<double>(get_J(), 0.0);
                        if ((i == k) & (j == l) & (i != j)) {
                            for (int p = 0; p < connections.at(i).size(); p++) {
                                if (j == connections.at(i).at(p)) {
                                    U_matrix[i][j][k][l] = complex<double>(get_J(), 0.0);
                                    /* Works */
                                    //cout << "Nonlocal: [" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                                }
                            }
                        }
                        if ((i == j) & (i == k) & (i == l)) {
                            U_matrix[i][j][k][l] = complex<double>(U, 0.0);
                            /* Works */
                            //cout << "Local: [" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                        }
                    }
                }
            }
        }
    }

    void print_U_matrix(){
        int lat_s = connections.size();
        cout << "connections size" << connections.size() << endl;
        cout << "U matrix:" << endl;
        for(int i = 0; i < lat_s; i++){
            for (int j = 0; j < lat_s; j++){
                for(int k = 0; k < lat_s; k++){
                    for(int l = 0; l < lat_s; l++){
                        if ((i == j) & (i == k) & (i == l)){
                            cout << "Local: \t\t[" << i << "][" << j << "][" << k << "][" << l << "] = " << U_matrix[i][j][k][l] << endl;
                        }
                        if ((i == k) & (j == l) & (i != j)){
                            for (int p = 0; p < connections.at(i).size(); p++) {
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

    complex<double> fermionic_matsubara(int n){
        return double(M_PI * (2 * n + 1) / get_beta()) * complex<double> (0, 1);
    }

    complex<double> bosonic_matsubara(int n){
        return double(M_PI * 2 * n / get_beta()) * complex<double> (0, 1);
    }

    double get_U() const {
        return U;
    }

    double get_J() const {
        return J;
    }

    double get_t() const {
        return t;
    }

    complex<double> get_chemical_potential() const {
        return complex<double> (mu, 0.0);
    }

    double get_beta() const {
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

    vector <vector <vector < vector <std::complex<double> > > > > get_U_matrix() const {
        return U_matrix;
    }

    vector < vector <std::complex < double > > > get_t_matrix() const {
        return t_matrix;
    }

    vector <complex<double>> get_omega() const {
        return omega;
    }

    vector <complex<double>> get_nu() const {
        return nu;
    }

private:

    void set_U(double local_coulomb){
        U = local_coulomb;
    }

    void set_J(double nonlocal_exchange){
        J = nonlocal_exchange;
    }

    void set_chemical_potential(double chemical_potential){
        mu = chemical_potential;
    }

    void set_hopping(double hopping){
        t = hopping;
    }

    void set_beta(double inversive_temperature){
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
