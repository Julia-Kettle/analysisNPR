#ifndef IO_H
#define IO_H

#include <iostream>
#include <sys/stat.h>
#include <string>
#include "Grid/Grid.h"
#include <linux/limits.h>
#include <sstream>

/////////////////////////////////////////////////////////////////////////////////////
// reads Nconf HDF5 files and returns vector of data
template<typename T>
void  readDataByConfig(std::string filestem, std::string groupLabel, std::vector<int> configs, std::vector<T> &data_vector)
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

/////////////////////////////////////////////////////////////////////////////
// overload Grids read function for int data types
void read(Grid::XmlReader &reader,std::string name, int &output)
{
    std::string tmp;
    read(reader,name,tmp);
    std::istringstream ss(tmp);
    ss >> output;
}


/////////////////////////////////////////////////////////////////////////////
// overload Grids read function for vector data types
template <typename T>
void read(Grid::XmlReader &reader,std::string name, std::vector<T> &output)
{
    std::string tmp;
    read(reader,name,tmp);
    std::stringstream iss( tmp );
    T number;
    output.resize(0);
    while( iss >> number )
    {
        output.push_back(number);
    }
}

////////////////////////////////////////////////////////////////////////////
// wraps the read function - returns template data type without having to provide
template <typename T>
T parseParam(Grid::XmlReader &reader,const std::string &label) 
{
    T data;
    read(reader, label, data);
    return data;
}


// recursive mkdir /////////////////////////////////////////////////////////////
//  coppied from Hadrons  - Antonin portelli's code. Should probably just import hadrons too
int mkdir(const std::string dirName)
{
    if (!dirName.empty() and access(dirName.c_str(), R_OK|W_OK|X_OK))
    {
        mode_t mode755;
        char   tmp[PATH_MAX];
        char   *p = NULL;
        size_t len;

        mode755 = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;

        snprintf(tmp, sizeof(tmp), "%s", dirName.c_str());
        len = strlen(tmp);
        if(tmp[len - 1] == '/')
        {
            tmp[len - 1] = 0;
        }
        for(p = tmp + 1; *p; p++)
        {
            if(*p == '/')
            {
                *p = 0;
                ::mkdir(tmp, mode755);
                *p = '/';
            }
        }

        return ::mkdir(tmp, mode755);
    }
    else
    {
        return 0;
    }
}

//wrapper for writing to file
template <typename T>
void save_result(std::string output_file, std::string label, T result)
{
    // Find the name of the path and create if doens't exist
    size_t index = output_file.find_last_of("/");
    std::string output_dir = output_file.substr(0,index);
    mkdir(output_dir);

    // Write hdf5 file
    Grid::Hdf5Writer writer(output_file);
    Grid::write(writer,label,result);
}

template <typename T>
void save_text(std::string output_file, std::string label, T result)
{
    // Find the name of the path and create if doens't exist
    size_t index = output_file.find_last_of("/");
    std::string output_dir = output_file.substr(0,index);
    mkdir(output_dir);

    // Write hdf5 file
    Grid::TextWriter writer(output_file);
    Grid::write(writer,label,result);

} 

#endif
