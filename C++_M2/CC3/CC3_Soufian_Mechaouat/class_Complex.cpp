#include "class_Complex.hpp"
#include <cmath>

const complex& complex::operator=(complex& z){
    x = z.X();
    y = z.Y();
    return *this;
}

const complex& complex::operator+=(complex& z ){
    x += z.X();
    y += z.Y();
    return *this;
}