//unveraendert von Thilo uebernommen
#ifndef MATRIX2_H
#define MATRIX2_H

#ifndef Matrixtyp
#define Matrixtyp double
#endif  /* Matrixtyp */

#include <stdio.h>
#include <iostream>  //.h weg
#include <iomanip>   //.h weg
#include <math.h>
#include <malloc.h>
#include "daten.h"
#include <string.h>

struct Vektor;

struct Matrix
{	
  void neue_Dimension(int Ze, int Sp, Matrixtyp Wert = 0);
  Matrixtyp& operator () (int Zeile, int Spalte);
  Matrix Teilmatrix(int oZ, int uZ, int lS, int rS);
  Matrix Zeilenmatrix(int Zeile)   { return Teilmatrix(Zeile,Zeile,1,Spalten);  };
  Matrix Spaltenmatrix(int Spalte) { return Teilmatrix(1,Zeilen,Spalte,Spalte); };
  Matrix Transponiert();
  Matrix operator = (const Matrix& M);
  Matrix operator = (Matrixtyp Wert);
  Matrix operator = (int Wert) { return operator =  ((Matrixtyp) Wert); };
  Matrix operator + (const Matrix& M);
  Matrix operator += (const Matrix& M);
  Matrix operator - (const Matrix& M);
  Matrix operator -= (const Matrix& M);
  Matrix operator * (Matrixtyp Wert);
  Matrix operator *= (Matrixtyp Argument);
  Matrix operator * (const Matrix& M);
  Matrix operator *= (const Matrix& M);
  Matrix operator / (Matrixtyp Wert);
  Matrix operator /= (Matrixtyp Nenner) { return operator *= (1/Nenner);};
  int Z() const { return Zeilen; };
  int S() const { return Spalten; };
  Matrixtyp Summe();
  Matrixtyp Norm();
  Matrixtyp** nr();
  Matrixtyp Max(int sig = 0);
  Matrixtyp Min(int sig = 0);

  friend std::ostream& operator << (std::ostream& s,const Matrix& M);
  friend std::istream& operator >> (std::istream& s,Matrix& M);

  Matrix(int Zeilen, int Spalten, Matrixtyp Wert = 0, int Z_O = 1, int S_O = 1);
  Matrix(int Zeilen, int Spalten, Matrixtyp **M, int Z_O = 1, int S_O = 1);
  Matrix(Matrixtyp **M, int Zeilen, int Spalten, int Z_O = 1, int S_O = 1);
  Matrix(const Matrix& M);
  ~Matrix();
  Matrixtyp **Elm;
  protected :
  int Zeilen, Spalten, Z_Offset, S_Offset;
  void Dimension_pruefen(const Matrix &M);
  Matrixtyp** Matrix_erzeugen(Matrixtyp Wert = 0, int Z = 0, int S = 0);
};

  std::ostream& operator << (std::ostream& s,const Matrix& M);
  std::istream& operator >> (std::istream& s,Matrix& M);

struct Vektor : Matrix
{	
  Vektor(int Elemente, Matrixtyp Wert = 0, int Offset = 1) :
    Matrix(1, Elemente, Wert, 1, Offset) {};
  Vektor(int Elemente, Matrixtyp **M, int Z_O = 1, int S_O = 1) :
	 Matrix(1, Elemente, M, Z_O, S_O) {};
  Vektor(int Elemente, Matrixtyp *Vektor) : Matrix(1, Elemente, &Vektor,0,1) {};
  Vektor(Matrixtyp *V, int Elemente, int Offset = 1) :
    Matrix(&V-1, 1, Elemente, 1, Offset) {};
  Vektor(Matrixtyp **V, int Elemente, int Offset = 1) :
    Matrix(V, 1, Elemente, 1, Offset) {};
  void neue_Dimension(int N, Matrixtyp Wert = 0) { Matrix::neue_Dimension(1, N, Wert); };
  Matrixtyp& operator () (int Element) { return Matrix::operator()(1,Element);};
  Matrixtyp operator* (Vektor& Multiplikant);
  Matrix operator* (Matrixtyp Multiplikant)
    { return Matrix::operator*(Multiplikant); };
  Matrix operator = (Matrixtyp Wert) { return Matrix::operator = (Wert); };
  Vektor operator = (const Matrix& M);
  Vektor& operator>>= (int rshift);
  Vektor& operator<<= (int lshift);
  Matrixtyp* nr() { return (Elm[0] - 1); };
  int E() { return Spalten; };
};


#if (Matrixtyp==double)
#ifndef _Mathe_to_h

  struct dm : Matrix
{
  dm(int Zeilen, int Spalten, Matrixtyp Wert = 0, int Z_O = 1, int S_O = 1) :
    Matrix(Zeilen, Spalten, Wert, Z_O, S_O) {};

};
  struct dv : Vektor
{	
  dv(int Elemente, Matrixtyp Wert = 0) : Vektor(Elemente, Wert) {};
};
/*
  struct edm : Matrix
  {
    edm(int n) : Matrix(n,n) { for (int i = 0; i < n; Elm[i][i] = 1, i++); };
  };
*/

#endif /* ndef _Mathe_to_h      */
#endif /* (Matrixtyp == double) */

#endif /*MATRIX_H*/

/* End. */
