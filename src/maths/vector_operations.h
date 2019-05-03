#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
    
//transpose vector of vectors
template<typename T>
std::vector<std::vector<T>> transpose(std::vector<std::vector<T>> vec_in)
{
    size_t n_rows = vec_in.size();
    size_t n_cols = vec_in[0].size();
    
    std::vector<std::vector<T>> vec_out(n_cols, std::vector<T>(n_rows));

    for(int i=0;i<n_rows;i++)
    for(int j=0;j<n_cols;j++)
    {
        vec_out[j][i] = vec_in[i][j];
    }
    return vec_out;
}

//Take real part of vector of complex
std::vector<Grid::Real> real(std::vector<Grid::ComplexD> comp_vec)
{
    std::vector<Grid::Real> reals;
    for( auto elem : comp_vec )
    {
        reals.push_back(real(elem));
    }
    return reals;
}


#endif
