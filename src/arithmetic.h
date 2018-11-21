#ifndef arithmetic_h
#define arithmetic_h

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

#endif
