#pragma once

#include <iostream>
#include "class_Point2D.hpp"

class complex: public point
{
  public:

   complex(double xx, double yy):point(xx,yy){};
   complex(const point&p ):point(p.X(),p.Y()){};

   const complex& operator=(complex& z);
   const complex& operator+=(complex& z );

};
