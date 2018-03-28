//
// Created by Cooper on 3/28/18.
//

#ifndef SECON_ORGER_DIAGRAM_SIGMA_GREENFUNCTION_H
#define SECON_ORGER_DIAGRAM_SIGMA_GREENFUNCTION_H

#include <complex>
#include "Data.h"
#include "boost/multi_array.hpp"
#include <cassert>

class GreenFunction {
private:
    int number_of_spins, lat_s, omega_length;

    typedef boost::multi_array<complex<double >, 4> array_type;
    typedef array_type::index index;
    array_type green_function;
public:
    void setGreen_function(const array_type &green_function) {
        GreenFunction::green_function = green_function;
    }

public:
    const array_type &getGreen_function() const {
        return green_function;
    }

public:
    GreenFunction(const Data &data){
        number_of_spins = data.get_number_of_spins();
        lat_s = data.get_number_of_sites();
        omega_length = data.get_number_of_freq();
        cout << "from GF class : " << number_of_spins << ", " << lat_s << ", " << omega_length << endl;
        array_type function(boost::extents[number_of_spins][lat_s][lat_s][omega_length]);
        for (int spin = 0; spin < number_of_spins; spin++) {
            for (int i = 0; i < lat_s; i++) {
                for (int j = 0; j < lat_s; j++) {
                    for (int freq = 0; freq < omega_length; freq++) {
                        function[spin][i][j][freq] = complex<double> (0.0, 0.0);
                    }
                }
            }
        }
        cout << function.shape();
        setGreen_function(function);
    }

};


#endif //SECON_ORGER_DIAGRAM_SIGMA_GREENFUNCTION_H
