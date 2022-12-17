#include <iostream>
#include "class_pointXd.hpp"

using namespace std;

int main()
{
  int    a[1]={1};
  float  b[2]={2., 2.};
  double c[3]={3., 3., 3.};
  double d[3]={4., 4., 4.};

  pointXd<int,   1> p1(a);
  pointXd<float, 2> p2(b);
  pointXd<double,3> p3(c);

  cout << "p1: " << endl;
  p1.print();
  cout << endl;

  cout << "p2: " << endl;
  p2.print();
  cout << endl;

  cout << "p3: " << endl;
  p3.print();
  cout << endl;

  // test for copy constructor
  pointXd<double,3> p4(p3);
  cout << "p4: " << endl;
  p4.print();
  cout << endl;

  // test for = constructor
  pointXd<double,3> p5(d);
  p5=p4;
  cout << "p5: " << endl;
  p5.print();
  cout << endl;
}
