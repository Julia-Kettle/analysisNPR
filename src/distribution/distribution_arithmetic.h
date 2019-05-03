#ifndef DISTRIBUTION_ARITHMETIC_H
#define DISTRIBUTION_ARITHMETIC_H


/////////////////////////////////// Operator Overloading ( left-wise ) /////////////////////////////////////////

// V op Distribution<T>

template<typename T, typename V>
Distribution<T> operator * (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs*rhs.get_values()); }

template<typename T, typename V>
Distribution<T> operator + (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs+rhs.get_values()); }

template<typename T, typename V>
Distribution<T> operator - (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs-rhs.get_values()); }

template<typename T, typename V>
Distribution<T> operator / (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs*(std::vector<T>(rhs.get_values().size(),1.0)/rhs.get_values())); }

#endif
