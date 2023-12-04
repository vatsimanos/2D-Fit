//unveraendert von Thilo uebernommen
#include "matrix2.h"

void Matrix::Dimension_pruefen(const Matrix& M)
{
  if ((Zeilen != M.Zeilen) || (Spalten != M.Spalten) ||
      (Z_Offset != M.Z_Offset) || (S_Offset != M.S_Offset))
  {
    std::cout << "Matrixoperation mit Matritzen unterschiedlicher Größe" << std::endl;
    exit(1);
  }
}

Matrixtyp** Matrix::Matrix_erzeugen(Matrixtyp Wert, int Z, int S)
{	
  int i, j;
  Matrixtyp** Neue_Matrix;
  if (Z == 0) Z = Zeilen;
  if (S == 0) S = Spalten;
 

  Neue_Matrix = (Matrixtyp**)malloc(sizeof(Matrixtyp*) * Z);
  for(i = 0; i < Z; i++)
  {
    Neue_Matrix[i] = (Matrixtyp*)malloc( sizeof(Matrixtyp) * S);
        for(j = 0; j < S; Neue_Matrix[i][j] = Wert, j++);
  }
  
  return Neue_Matrix;
}

Matrix Matrix::Teilmatrix(int oZ, int uZ, int lS, int rS)
{
  Matrixtyp **Neue_Matrix;
  int i, j, Ze, Sp;

  if (oZ < Z_Offset || uZ >= Zeilen + Z_Offset ||
      lS < S_Offset || rS >= Spalten + S_Offset)
  {
    std::cout << "Unzulässige Indizes bei Funktion : Teilmatrix()" << std::endl;
    exit(1);
  }
  Ze = uZ - oZ + 1;
  Sp = rS - lS + 1;
  oZ -= Z_Offset;
  lS -= S_Offset;
  Neue_Matrix = Matrix_erzeugen(0, Ze, Sp);
  for(i = 0; i < Ze; i++)
   for(j = 0; j < Sp; j++)
    Neue_Matrix[i][j] = Elm[i + oZ][j + lS];
  return Matrix(Ze, Sp, Neue_Matrix);
}

Matrix Matrix::Transponiert()
{
  int i, j;
  Matrixtyp **Neue_Matrix;
  Neue_Matrix = Matrix_erzeugen();
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Neue_Matrix[i][j] = Elm[j][i], j++);
  return Matrix(Zeilen, Spalten, Neue_Matrix);
}

void Matrix::neue_Dimension(int Ze, int Sp, Matrixtyp Wert)
{
  int i;
  for(i = 0; i < Zeilen; free((char *)Elm[i]), i++);
  free((char *)Elm);
  Zeilen = Ze;
  Spalten = Sp;
  Elm = Matrix_erzeugen(Wert);
}

Matrix Matrix::operator = (const Matrix &M)
{
  int i, j;
  Dimension_pruefen(M);
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; j++)
    Elm[i][j] = M.Elm[i][j];
  return *this;
}

Matrix Matrix::operator = (Matrixtyp Wert)
{
  int i, j;
  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; Elm[i][j] = Wert, j++);
  return *this;
}

Matrixtyp& Matrix::operator () (int Zeile, int Spalte)
{
  if (Zeile >= Z_Offset && Zeile < Zeilen + Z_Offset
      && Spalte >= S_Offset && Spalte < Spalten + S_Offset)
    return Elm[Zeile - Z_Offset][Spalte - S_Offset];
  else { std::cerr << "Zugriff auf nicht vorhandenes Matrixelement\n"; abort(); }
}

Matrix Matrix::operator + (const Matrix& M)
{
  int i, j;
  Matrixtyp **Neue_Matrix;

  Dimension_pruefen(M);
  Neue_Matrix = Matrix_erzeugen();
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Neue_Matrix[i][j] = Elm[i][j] + M.Elm[i][j], j++);
  return Matrix(Zeilen, Spalten, Neue_Matrix);
}

Matrix Matrix::operator += (const Matrix& M)
{
  int i, j;
  Dimension_pruefen(M);
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Elm[i][j] = Elm[i][j] + M.Elm[i][j], j++);
  return *this;
}

Matrix Matrix::operator - (const Matrix& M)
{
  int i, j;
  Matrixtyp** Neue_Matrix;

  Dimension_pruefen(M);
  Neue_Matrix = Matrix_erzeugen();
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Neue_Matrix[i][j] = Elm[i][j] - M.Elm[i][j], j++);
  return Matrix(Zeilen, Spalten, Neue_Matrix);
}

Matrix Matrix::operator -= (const Matrix& M)
{
  int i, j;
  Dimension_pruefen(M);
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Elm[i][j] = Elm[i][j] - M.Elm[i][j], j++);
  return *this;
}

Matrix Matrix::operator * (Matrixtyp Argument)
{
  int i, j;
  Matrixtyp** Neue_Matrix;
  Neue_Matrix = Matrix_erzeugen();
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Neue_Matrix[i][j] = Elm[i][j] * Argument, j++);
  return Matrix(Zeilen, Spalten, Neue_Matrix);
}

Matrix Matrix::operator * (const Matrix& M)
{
  int i, j, k, transponiert, neu_Spalten;
  Matrixtyp** Neue_Matrix;

  if ((transponiert = (Spalten == M.Spalten)) || (Spalten == M.Zeilen))
  {
    if (transponiert) neu_Spalten = M.Zeilen;
      else neu_Spalten = M.Spalten;    
    Neue_Matrix = Matrix_erzeugen(0,Zeilen,neu_Spalten);
    for(i = 0; i < neu_Spalten; i++)
     for(j = 0; j < Zeilen; j++)
      for(k = 0; k < Spalten; k++) Neue_Matrix[j][i] +=
        Elm[j][k] * (transponiert ? M.Elm[i][k] : M.Elm[k][i]);
  }
  else
  {
    std::cout<<"Matrixmultiplikation mit falschen Dimensionen"<<std::endl;
    exit(1);
  }
  return Matrix(Zeilen, neu_Spalten, Neue_Matrix);
}

Matrix Matrix::operator *= (Matrixtyp Argument)
{
  int i, j;
  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; j++) Elm[i][j] *= Argument;
  return *this;
}

Matrix Matrix::operator / (Matrixtyp Nenner)
{
  int i, j;
  Matrixtyp** Neue_Matrix;
  Neue_Matrix = Matrix_erzeugen();
  for(i = 0; i < Zeilen; i++)
   for(j = 0; j < Spalten; Neue_Matrix[i][j] = Elm[i][j] / Nenner, j++);
  return Matrix(Zeilen, Spalten, Neue_Matrix);
}

Matrixtyp Matrix::Summe()
{
  int i,j;
  double Sum = 0;

  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; j++) Sum += Elm[i][j];
  return Sum;
}

Matrixtyp Matrix::Norm()
{
  int i,j;
  double Norm = 0;

  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; j++) Norm += (Elm[i][j] * Elm[i][j]);
  return sqrt(Norm);
}

Matrixtyp** Matrix::nr()
{
  int i;
  Matrixtyp **nrM;
  nrM = (Matrixtyp **)malloc(Zeilen * sizeof(Matrixtyp *));
  nrM -= Z_Offset;
  for(i = 0; i < Zeilen; nrM[i + Z_Offset] = Elm[i] - S_Offset, i++);
  return nrM;
}

Matrixtyp Matrix::Max(int sig)
{
  int i, j, Flag = 1;
  Matrixtyp Wert = 0, Zw;

  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; j++)
    {
      if ((((Zw = Elm[i][j] * sig) > 0) || (sig == 0)) && (1 == Flag--))
        Wert = Zw;
      if ((Wert <= (Zw = Elm[i][j])) && ((Zw * sig > 0) || (sig == 0)))
        Wert = Zw;
    }
  return Wert;
}

Matrixtyp Matrix::Min(int sig)
{
  int i, j, Flag = 1;
  Matrixtyp Wert = 0, Zw;

  for (i = 0; i < Zeilen; i++)
    for (j = 0; j < Spalten; j++)
    {
      if ((((Zw = Elm[i][j] * sig) > 0) || (sig == 0)) && (1 == Flag--))
        Wert = Zw;
      if ((Wert >= (Zw = Elm[i][j])) && ((Zw * sig > 0) || (sig == 0)))
        Wert = Zw;
    }
  return Wert;
}

Matrix::Matrix(int Z, int S, Matrixtyp Wert, int Z_O, int S_O)
{		
  Zeilen = Z;
  Spalten = S;
  Z_Offset = Z_O;
  S_Offset = S_O;
  Elm = Matrix_erzeugen(Wert);
}

Matrix::Matrix(int Z, int S, Matrixtyp** M, int Z_O, int S_O)
{  			
  Zeilen = Z;
  Spalten = S;
  Z_Offset = Z_O;
  S_Offset = S_O;
  Elm = M;
}

Matrix::Matrix(Matrixtyp **M, int Z, int S, int Z_O, int S_O)
{			
  int i, j;
  Zeilen = Z;
  Spalten = S;
  Z_Offset = Z_O;
  S_Offset = S_O;
  Elm = Matrix_erzeugen();
    for(i = 0; i < Zeilen; i++)
      for(j = 0; j < Spalten; j++) Elm[i][j] = M[i+Z_O][j+S_O];
}

Matrix::Matrix(const Matrix& M)
{		
  int i, j;
  Zeilen = M.Zeilen;
  Spalten = M.Spalten;
  Z_Offset = M.Z_Offset;
  S_Offset = M.S_Offset;
  Elm = Matrix_erzeugen();
    for(i = 0; i < Zeilen; i++)
      for(j = 0; j < Spalten; Elm[i][j] = M.Elm[i][j], j++);
}


Matrix::~Matrix()
{	
  int i;
  for(i = 0; i < Zeilen; i++)
    if (Elm[i]!=NULL)
	free((char *)Elm[i]);
        
  if (*Elm!=NULL)
       free((char *)Elm);        

}

std::ostream& operator << (std::ostream& s, const Matrix& M)
{
  int i, j;

  for (i = 0; i < M.Zeilen; std::cout << std::endl, i++)
   for (j = 0; j < M.Spalten; s << "  " << M.Elm[i][j], j++);
  return s;
}

std::istream& operator >> (std::istream& s, Matrix& M)
{
  int i, j, Z, S;

  s >> Z >> S;
  M.neue_Dimension(Z,S);
  for(i = 1; i <= Z; i++)
    for(j = 1; j <= S; s >> M(i,j), j++);
  return s;
}

Vektor Vektor::operator = (const Matrix& M)
{
  int i, transponiert;	
  if ((transponiert = (M.Z() == Spalten && M.S() == 1)) ||
       (M.Z() == 1 && M.S() == Spalten))
    for (i = 0; i < Spalten; i++)
      Elm[0][i] = transponiert ? M.Elm[i][0] : M.Elm[0][i];
    else { 
    	  std::cout << "Falsche Dimension bei Vektormanipulation\n"; 
    	  exit(1); 
    	 }
  return *this;
}

Vektor& Vektor::operator >>= (int rshift)
{
  int i;
  for (i=Spalten-1; i>=rshift; i--) Elm[0][i] = Elm[0][i-rshift];
  for (i=0; i<rshift; i++) Elm[0][i] = 0;
  return *this;
}

Vektor& Vektor::operator <<= (int lshift)
{
  int i;
  for (i=0; i<Spalten-lshift; i++) Elm[0][i] = Elm[0][i+lshift];
  for (i=Spalten-lshift; i<Spalten; i++) Elm[0][i] = 0;
  return *this;
}

Matrixtyp Vektor::operator* (Vektor& Multiplikant)
{
  int i;
  Matrixtyp dummy = 0;
  Dimension_pruefen(Multiplikant);
  for (i=0; i<Spalten; i++) dummy += Elm[0][i] * Multiplikant.Elm[0][i];
  return dummy;
}


