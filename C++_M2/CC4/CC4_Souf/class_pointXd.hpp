#pragma once
#include<iostream>

using namespace std;

template <typename S, int Nd>
class pointXd{

protected:

  S* coord;

public:

  pointXd(S* vv);
  inline S x() const{return coord[0];};
  inline S y() const{return coord[1];};
  void print() const;
  ~pointXd();  //  destructor
  pointXd(const pointXd& p);  //  copy constructor
  pointXd<S,Nd> & operator =(pointXd<S,Nd> const &vv);
};

//Def

template <typename S, int Nd>
void pointXd<S,Nd>::print() const{
  if(Nd==1){cout << "x = " << coord[0]<<endl;}
  else if(Nd==2){cout << "(x,y) = (" << coord[0]<<","<<coord[1]<<")"<<endl;}
  else if(Nd==3){cout << "(x,y,z) = (" << coord[0]<<","<<coord[1]<<","<<coord[2]<<")"<<endl;}
};


template <typename S, int Nd>
pointXd<S,Nd>::pointXd(S* vv){
  coord = new S[Nd];
  for (int i=0; i<Nd; i++){
    coord[i]=vv[i];
  }
};

template <typename S, int Nd>
pointXd<S,Nd>::~pointXd(){delete[] coord;};

template <typename S, int Nd>
pointXd<S,Nd>::pointXd(const pointXd& p){
  coord = new S[Nd];
  for (int i=0; i<Nd; i++){
    coord[i]=p.coord[i];
  }
}  

//3
template<typename S, int Nd>
pointXd<S,Nd>& pointXd<S,Nd>::operator =(pointXd<S,Nd> const & vv){
    if(this != &vv){
        delete [] coord;
        coord = new S[Nd];
        for(int i=0; i<Nd; ++i){
            coord[i]=vv.coord[i];
        }
    }
    return *this;
}