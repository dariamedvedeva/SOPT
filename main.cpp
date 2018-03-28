#include <iostream>

#include "Data.h"
#include "Calculation.h"

using namespace std;

int main() {
    cout << "ELLO MY QUEEN!" << endl;

    int iteration = 0;
    bool half_fill;
    half_fill = true;
    /* Parameters */
    double local_coulomb         = 1.0;
    double nonlocal_exchange     = 0.05;
    double hopping               = -0.25;
    double inversive_temperature = 10.0;
    int number_of_frequencies   = 10; // number of frequencies
    int number_of_sites         = 4;

    double chemical_potential;

    if (half_fill){
        chemical_potential    = 0.0;
    } else {
        chemical_potential    = local_coulomb / 2.0;
    }

    Data data;
    data.init_parameters(local_coulomb, nonlocal_exchange,
                         chemical_potential, hopping,
                         inversive_temperature, number_of_frequencies,
                         number_of_sites);
    data.print_data();

    data.frequencies();
    data.construct_connections();
    data.print_connections();

    data.construct_hopping_matrix();
    data.print_hopping_matrix();

    data.construct_U_matrix();
    data.print_U_matrix();

    /* Initial lattice Green's function */
    Calculation calculation(data);

    cout << "Zero iteration." << endl;

    /* G0 */
    calculation.construct_initial_lattice_function();

    /* GF = G0 */
    calculation.GF_takes_GO();

    /* Sigma */
    cout << "Sigma computing..." << endl;
    calculation.compute_Sigma();

    /* G_final */
    cout << "Dyson equation..." << endl;
    calculation.Dyson_equation();

    cout << "Start cycle" << endl;
    while(calculation.test_convergency() == false){
        iteration += 1;
        cout << "Iteration number " << iteration << endl;

        /* GF = GF_final*/
        calculation.GF_takes_GF_final(); /* and inversive matrix too*/

        /* Sigma */
        cout << "Sigma computing..." << endl;
        calculation.compute_Sigma();
        /* GF_final */
        cout << "Dyson equation..." << endl;
        calculation.Dyson_equation();
    }


    calculation.print_Sigma_in_file();
    calculation.print_GF_in_file();
    calculation.print_GF0_in_file();
    cout << "Total number of iterations = " << iteration << endl;
    data.clear_memory_t_matrix();

    /*Memory clean block*/
    //calculation.clean_memory_GF();
    cout << "END" << endl;
    return 0;
}

