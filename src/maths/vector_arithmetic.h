#ifndef VECTOR_ARITHMETIC_h
#define VECTOR_ARITHMETIC_h

#include "distribution.h"
#include "Grid/Grid.h"

// c = a op b -> c[i] = a[i] op b[i[
// op = {+,-,/,*}

///////////Binary operators////////////////   
// vector + vector  
template <typename T>
std::vector<T> operator +(std::vector<T> a, std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]+b[i]); }
    return c;
}

template <typename T>
std::vector<T> operator -( std::vector<T> a,  std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]-b[i]); }
    return c;
}

template <typename T>
std::vector<T> operator *( std::vector<T> a,  std::vector<T> b)
{
    //if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]*b[i]); }
    return c;
}

template <typename T>
std::vector<T> operator /( std::vector<T> a,  std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]/b[i]); }
    return c;
}


///////////Binary operators////////////////   
// vector + vector  
template <typename T, typename U>
std::vector<T> operator +(std::vector<U> a, std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]+b[i]); }
    return c;
}

template <typename T, typename U>
std::vector<T> operator -( std::vector<U> a,  std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]-b[i]); }
    return c;
}

template <typename T, typename U>
std::vector<U> operator *( std::vector<U> a,  std::vector<T> b)
{
    //if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]*b[i]); }
    return c;
}

template <typename T, typename U>
std::vector<T> operator /( std::vector<U> a,  std::vector<T> b)
{
    if(a.size() != b.size()){ throw std::string("vectors must be of equal size"); }
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]/b[i]); }
    return c;
}




///////////Binary operators////////////////   
// vector + scalar 
template <typename T, typename S>
std::vector<T> operator +( std::vector<T> a,  S b)
{
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]+b); }
    return c;
}

template <typename T, typename S>
std::vector<T> operator -( std::vector<T> a,  S b)
{
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]-b); }
    return c;
}

template <typename T, typename S>
std::vector<T> operator *( std::vector<T> a,  S b)
{
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]*b); }
    return c;
}

    template <typename T, typename S>
std::vector<T> operator *(S a, std::vector<T> b)
{
    std::vector<T> c;
    c.reserve(b.size());
    for (int i=0; i< b.size(); i++){ c.push_back(a*b[i]); }
    return c;
}

template <typename T, typename S>
std::vector<T> operator /( std::vector<T> a,  S b)
{
    std::vector<T> c;
    c.reserve(a.size());
    for (int i=0; i< a.size(); i++){ c.push_back(a[i]/b); }
    return c;
}

//sqrt
template <typename T>
std::vector<T> sqrt(std::vector<T> input)
{
    std::vector<T> result;
    for (auto elem : input)
    {
        auto tmp = sqrt(elem);
        result.push_back(tmp);
    }
    return result;
}

#endif
