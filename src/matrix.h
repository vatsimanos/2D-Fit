/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *
 *  templateklassen matrix, tensor als aufsatz zum stl_vector
 *  nur speicherverwaltung, keine  speziellen operatoren!
 *  von Oliver Radomski, August 1999
 * 
 *  bei aufruf der konstruktoren matrix<>(), tensor<>()
 *  wird noch kein speicher bereitgestellt! resize(...) noetig!
 *  oder matrix<>(x,y)
 *******************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector> //.h weg

template <class T> class matrix
{
  std::vector< std::vector<T> > mat;   //std:: added und folgende
public:
  matrix() {}
  matrix(size_t x, size_t y) {resize(x, y, T());}
  matrix(size_t x, size_t y, const T& t) {resize(x, y, t);}
  ~matrix() {}

  // mat[x][y]
  size_t sizex() const {return mat.size();}
  size_t sizey() const {return (mat.size()) ? mat[0].size() : 0 ;}

  std::vector<T>& operator[](size_t i) {return mat[i];}
  const std::vector<T>& operator[](size_t i) const {return mat[i];}

  void resize(size_t x, size_t y, const T& t);
  void resize(size_t x, size_t y) {resize(x, y, T());}
  void clear();

  bool empty() const {return mat.empty();}
};

template <class T> class tensor
{
  std::vector< matrix<T> > ten;
public:
  tensor() {}
  tensor(size_t x, size_t y, size_t z) {resize(x, y, z, T());}
  tensor(size_t x, size_t y, size_t z, const T& t = T()) {resize(x, y, z, t);}
  ~tensor() {}

  // ten[x][y][z]
  size_t sizex() const {return ten.size();}
  size_t sizey() const {return (ten.size()) ? ten[0].sizex() : 0;}
  size_t sizez() const {return (ten.size()) ? ten[0].sizey() : 0;}

  matrix<T>& operator[](size_t i) {return ten[i];}
  const matrix<T>& operator[](size_t i) const {return ten[i];}

  void resize(size_t x, size_t y, size_t z, const T& t = T());
  void resize(size_t x, size_t y, size_t z) {resize(x, y, z, T());}
  void clear();

  bool empty() const {return ten.empty();}
};


template <class T> inline void 
matrix<T>::resize(size_t x, size_t y, const T& t)	
{
  if (x < sizex())
  {
    for (size_t i=0;i < x;i++) mat[i].resize(y, t);
    for (size_t i=x;i < sizex();i++) mat[i].clear();
    mat.resize(x);
  } else {
    mat.resize(x);
    for (size_t i=0;i < x;i++) mat[i].resize(y, t);
  }
}

template <class T> inline void 
matrix<T>::clear()
{
  for (size_t i=0;i < sizex();i++) mat[i].clear();
  mat.clear();
}

template <class T> void 
tensor<T>::resize(size_t x, size_t y, size_t z, const T& t)	
{
  if (x < sizex())
  {
    for (size_t i=0;i < x;i++) ten[i].resize(y, z, t);
    for (size_t i=x;i < sizex();i++) ten[i].clear();
    ten.resize(x);
  } else {
    ten.resize(x);
    for (size_t i=0;i < x;i++) ten[i].resize(y, z, t);
  }
}

template <class T> void 
tensor<T>::clear()
{
  for (size_t i=0;i < sizex();i++) ten[i].clear();
  ten.clear();
}

#endif /* _MATRIX_H_ */
