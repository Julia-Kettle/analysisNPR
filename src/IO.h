#ifndef IO_H
#define IO_H

#include <iostream>
#include <string>
#include "Grid/Grid.h"

/* This function reads in data for several configs 
 * Reads from HDF5 file format
 * results go into a data_vector of unknown type
 */

// read into vector
template<typename T>
void  readConfigs(std::string filestem, std::string groupLabel, std::vector<int> configs, std::vector<T> &data_vector)
{   
    
    data_vector.resize(configs.size()); // ensure data vector length correct
    
    // Loop through configs for reading 
    for(int iconf = 0; iconf<configs.size(); iconf++)
    {
        std::string filename = filestem + "." + std::to_string(configs[iconf]) +".h5";
        Grid::Hdf5Reader reader(filename);
        read(reader ,groupLabel, data_vector[iconf]);
    }
}




#endif