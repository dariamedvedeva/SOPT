#include <iostream>

#include "Data.h"
#include "Calculation.h"

using namespace std;

int main() {
    cout << "ELLO MY QUEEN!" << endl;

    int iteration = 0;

    /* Parametrs */
    float local_coulomb         = 1.0;
    float nonlocal_exchange     = 0.05;
    float chemical_potential    = local_coulomb / 2.0;
    float hopping               = -0.25;
    float inversive_temperature = 10.0;
    int number_of_frequencies   = 10; // number of frequencies
    int number_of_sites         = 4;

    Data data;
    data.init_parameters(local_coulomb, nonlocal_exchange,
                         chemical_potential, hopping,
                         inversive_temperature, number_of_frequencies,
                         number_of_sites);
    data.print_data();

    data.frequencies();
    data.construct_connections();

    data.construct_hopping_matrix();
    data.print_hopping_matrix();

    data.construct_U_matrix();
    data.print_U_matrix();

    /* Initial lattice Green's function */
    Calculation calculation(data);

    cout << "Zero iteration." << endl;

    /* G0 */
    calculation.construct_initial_lattice_function();
    calculation.inverse_matrix_to_G0();

    /* GF = G0 */
    calculation.GF_takes_GO();
    calculation.inverse_matrix_GF();

    /* Sigma */
    cout << "Sigma computing..." << endl;
    calculation.compute_Sigma();

    /* G_final */
    cout << "Dyson equation..." << endl;
    calculation.Dyson_equation();
    calculation.inverse_matrix_to_GF_final();

    cout << "Start cycle" << endl;
    int number_of_iterations = 10;

    while(calculation.test_convergency() == false){
        cout << "Iteration number " << iteration << endl;

        /* GF = GF_final*/
        calculation.GF_takes_GF_final(); /* and inversive matrix too*/

        /* Sigma */
        cout << "Sigma computing..." << endl;
        calculation.compute_Sigma();

        /* GF_final */
        cout << "Dyson equation..." << endl;
        calculation.Dyson_equation();
        calculation.inverse_matrix_to_GF_final();

        calculation.print_Sigma_in_file();
        calculation.print_GF_in_file();
        iteration += 1;

    }

    calculation.print_GF0_in_file();
    cout << "Total number of iterations = " << iteration << endl;

    /*Memory clean block*/
    data.clear_memory_t_matrix();
    //calculation.clean_memory_GF();

    cout << "END" << endl;
    return 0;
}

