/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <math.h>
#include "fit.h"
#include "gtkapp.h"
#include "error.h"
#include "y_shut.h"
#include "steady.h"

// init
Tlfit::Tlfit()
{
short i;

  singularMatrix=false;
  timeconstantSero=false;
  NumPoints=0;
  for (i=1;i < fit_dim+1;i++) ap[i]=0;
  for (i=1;i < map+1;i++) aa[i]=0;
}

// Least Square Routine
void 
Tlfit::lfit(RealArrayNData    x, 
            RealArrayNData    y,
            RealArrayNData    sig,             
            short            ndata, 
            RealArrayMA       a,                           
            short            ma,
            IntegerArrayMFIT  lista,
            short            mfit,
            RealArrayMAbyMA   covar,
            double           &chisq)
{
short           k,kk,j,ihit,i;
double          ym,wt,sum,sig2i;
RealArrayNPbyMP  beta;
RealArrayMA      afunc;
		 

  // array mit Null initialiesieren?
  kk = mfit+1;

  for(j=1;j < ma+1;j++)
  {
    ihit = 0;
    for (k=1;k < mfit+1;k++)
    {
      if (lista[k] == j) ihit = ihit+1;
    }
    if (ihit == 0)
    {
      lista[kk] = j;
      kk = kk+1;
    } else {
      if (ihit > 1)
      {
        warning((char *)"pause in routine LFIT\nimproper permutation in LISTA");
      }
    }
  }

  if (kk != ma+1)
  {
    warning((char *)"pause in routine LFIT\nimproper permutation in LISTA");
  }
  for (j=1;j < mfit+1;j++)
  {
    for (k=1;k < mfit+1;k++) covar[j][k] = 0.0;
    beta[j][1] = 0.0;
  }

  for (i=1;i < ndata+1;i++)
  {
    funcs(x[i],afunc);
    if (timeconstantSero) return;
    ym = y[i];
    if (mfit < ma)
    {
      for (j=mfit+1;j < ma+1;j++) ym = ym - a[lista[j]] * afunc[lista[j]];
      
    }
    sig2i = 1.0 / (sig[i] * sig[i]);
    for (j=1;j < mfit+1;j++)
    {
      wt = afunc[lista[j]] * sig2i;
      for (k=1;k < j+1;k++) covar[j][k] = covar[j][k] + wt * afunc[lista[k]];
      beta[j][1] = beta[j][1] + ym * wt;
    }
  }

  if (mfit > 1)
  {
    for (j= 2;j <  mfit+1;j++)
    {
      for (k=1;k < j;k++) covar[k][j] = covar[j][k];
    }
  }

  gaussj(covar,mfit,beta,1);

  if (singularMatrix) return;
  
  // keine negativen amplituden fuer gauss verteilungen
  for (j=1;j < mfit+1;j++) a[lista[j]] = ((beta[j][1] > 0) ? beta[j][1] : 0);

  chisq = 0.0;
  for (i=1;i < ndata+1;i++)
  {
    funcs(x[i],afunc);
    if (timeconstantSero) return;
    sum = 0.0;
    for (j=1;j <  ma+1;j++)
	{
	  sum += a[j] * afunc[j];
	}
    chisq += ((y[i] - sum) / sig[i]) * ((y[i] - sum) / sig[i]);
  }
  covsrt(covar,ma,lista,mfit);
}

double 
Tlfit::func(glnp pr)
{
// Fehlerfunkion,
// Quadratische Abweichung der Gaussumme von den Messdaten

short  i;

for (i=1;i < fit_dim+1;i++) ap[i] = pr[i];

  lfit(ax, ay, asig, NumPoints, aa, ama, alista, amfit, acovar, achisq);

  return achisq;
}

void
Tlfit::funcs(double z, RealArrayMA afunc)
{
// nur dummy fkt
}

void 
Tlfit::gaussj(RealArrayNPbyNP a,   
              short          n,
              RealArrayNPbyMP b,      
              short          m)
{
double         big=0,dum=0,pivinv=0;
short          icol=0,irow=0;
short          i,j,k,l,ll;
IntegerArrayNP  indxc;
IntegerArrayNP  indxr;
IntegerArrayNP  ipiv;           // IntegerArrayNP

  for (j=1;j <= n;j++) ipiv[j] = 0;
  for (i=1;i <= n;i++)
  {
    big = 0.0;
    for (j=1;j <= n;j++)
    {
      if (ipiv[j] != 1)
      {
        for (k=1;k <= n;k++)
        {
          if (ipiv[k] == 0)
          {
            if (abs(a[j][k]) >= big)
            {
              big = abs(a[j][k]);
              irow = j;
              icol = k;
            } else {
              if (ipiv[k] > 1)
              {
                if (!(Application.autom_steuer)) 
                {
                  warning((char *)"pause 1 in GAUSSJ\nsingular matrix");
                }
                singularMatrix = true;
                return;
              }
            }
          }
        }
      }  
    }
    ipiv[icol] = ipiv[icol]+1;
    if (irow != icol)
    {
      for (l=1;l <= n;l++)
      {
        dum = a[irow][l];
        a[irow][l] = a[icol][l];
        a[icol][l] = dum;
      }
      for (l=1;l <= m;l++)
      {
        dum = b[irow][l];
        b[irow][l] = b[icol][l];
        b[icol][l] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0)
    {
      if (!(Application.autom_steuer))
      {
        warning((char *)"pause 2 in GAUSSJ\nsingular matrix");
      } 
      singularMatrix = true;
      return;
    }
    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for (l=1;l <= n;l++) a[icol][l] = a[icol][l]*pivinv;
    for (l=1;l <= m;l++) b[icol][l] = b[icol][l]*pivinv;
    for (ll=1;ll <= n;ll++)
    {
      if (ll != icol) 
      {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l=1;l <= n;l++) a[ll][l] = a[ll][l]-a[icol][l]*dum;
        for (l=1;l <= m;l++) b[ll][l] = b[ll][l]-b[icol][l]*dum;
      }
    }
  }
  for (l=n;l >= 1;l--)
  {
    if (indxr[l] != indxc[l])
    {
      for (k=1;k <= n;k++)
      {
        dum = a[k][indxr[l]];
        a[k][indxr[l]] = a[k][indxc[l]];
        a[k][indxc[l]] = dum;
      }
    }
  }
}

void 
Tlfit::covsrt(RealArrayMAbyMA  covar,
              short           ma,
              IntegerArrayMFIT lista,
              short           mfit)
{
short  j,i;
double swap;

  for (j=1;j < ma;j++)
  {
    for (i=j+1;i < ma+1;i++) covar[i][j] = 0.0;
  }

  for (i=1;i < mfit;i++)
  {
    for (j=i+1;j < mfit+1;j++)
    {
      if (lista[j] > lista[i])
      {
        covar[lista[j]][lista[i]] = covar[i][j];
      } else {
        covar[lista[i]][lista[j]] = covar[i][j];
      }
    }
  }

  swap = covar[1][1];
  for (j=1;j < ma+1;j++)
  {
    covar[1][j] = covar[j][j];
    covar[j][j] = 0.0;
  }
  covar[lista[1]][lista[1]] = swap;
  for (j=2;j < mfit+1;j++) covar[lista[j]][lista[j]] = covar[1][j];
  for (j=2;j < ma+1;j++)
  {
    for (i=1;i < j;i++) covar[i][j] = covar[j][i];
  }
}

// Speicher belegen
Tamoeba::Tamoeba()
{
  ende = false;
}

// Downhill Simplex-Routine
void 
Tamoeba::amoeba(glmpnp  p,   
                glmp    y,
                short  ndim,
                double ftol,
                short  &iter)
{
const double alpha=1.0;
const double beta=0.5;
const double gamma=2.0;
const double itmax=1000;

short  mpts,j,inhi,ilo,ihi,i;
double yprr,ypr,rtol;
glnp    pr;
glnp    prr;
glnp    pbar;

  mpts = ndim+1;
  iter = 0;
  while (true)
  {
    ilo = 1;
    if (y[1] > y[2])
    {
      ihi = 1;
      inhi = 2;
    } else { 
      ihi = 2;
      inhi = 1;
    }
    for (i=1;i <= mpts;i++)
    {
      if (y[i] < y[ilo])  ilo = i;
      if (y[i] > y[ihi])
      {
        inhi = ihi;
        ihi = i;
      } else {
        if (y[i] > y[inhi])
        {
          if (i != ihi)  inhi = i;
        }
      }
    }
    rtol = 2.0 * fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) return;
    if (iter == itmax)
    {
      if (!(Application.autom_steuer))
      {
        warning((char *)"pause in AMOEBA\ntoo many iterations");
      }
      ende=true;
      return;
    }

//    if (!(Application.autom_steuer))
//    {
//      if (PeekMessage (Msg, HWindow, WM_RButtonDown, WM_RButtonDown, PM_Remove))
//      {
//        ende=true;
//        return;
//      }
//      SendMessage(HWindow, WM_Paint, 0,0);
//    }

    iter++;
    for (j=1;j <= ndim;j++) pbar[j] = 0.0;
    for (i=1;i <= mpts;i++)
    {
      if (i != ihi) for (j=1;j <= ndim;j++) pbar[j] = pbar[j]+p[i][j];
    }
    for (j=1;j <= ndim;j++)
    {
      pbar[j] = pbar[j]/ndim;
      pr[j] = (1.0+alpha)*pbar[j]-alpha*p[ihi][j];
    }
    ypr = func(pr);
    if ((singularMatrix) || (timeconstantSero)) return;

    if (ypr <= y[ilo])
    {
      for (j=1;j <= ndim;j++) prr[j] = gamma*pr[j]+(1.0-gamma)*pbar[j];
      yprr = func(prr);
      if ((singularMatrix) || (timeconstantSero)) return;
      if (yprr < y[ilo])
      {
        for (j=1;j <= ndim;j++) p[ihi][j] = prr[j];
        y[ihi] = yprr;
      } else {
        for  (j=1;j <= ndim;j++)  p[ihi][j] = pr[j];
        y[ihi] = ypr;
      }
    } else {
      if  (ypr >= y[inhi])
      {
        if (ypr < y[ihi])
        {
          for (j=1;j <= ndim;j++) p[ihi][j] = pr[j];
          y[ihi] = ypr;
        }
        for(j=1;j <= ndim;j++) prr[j] = beta*p[ihi][j]+(1.0-beta)*pbar[j];
        yprr = func(prr);
        if ((singularMatrix) || (timeconstantSero)) return;
        if (yprr < y[ihi])
        {
          for (j=1;j <= ndim;j++) p[ihi][j] = prr[j];
          y[ihi] = yprr;
        } else {
          for (i=1;i <= mpts;i++)
          {
            if (i != ilo)
            {
              for (j=1;j <= ndim;j++)
              {
                pr[j] = 0.5*(p[i][j]+p[ilo][j]);
                p[i][j] = pr[j];
              }
              y[i] = func(pr);
              if ((singularMatrix) || (timeconstantSero)) return;
            }
          }
        }
      } else {
        for (j=1;j <= ndim;j++) p[ihi][j] = pr[j];
        y[ihi] = ypr;
      }
    }
  }
}

// kombinierter Fit: Gauss-> Amplituden Simplex-> Kanalstr�me u. Abweichung
// double Tamoeba::fit(...)

void 
Tamoeba::startwerte()
{
short  i,j;

  for (i=1;i <= bndim;i++)
  {
    for (j=1;j <= bndim+1;j++)
	{
	  bp[j][i] = ap[i]; // Startwerte berechnen
	}
  }
  for (i=1;i <= bndim;i++)
  {
    bp[i][i] = bp[i][i] + 0.1 * ap[i]; // noch etwas addieren, damit nicht alle gleich
  }
}

void 
Tamoeba::funktionswerte(glnp cp)
{
int  i,j;

  for (i=1;i <= bndim+1;i++)
  {
    for (j=1;j <= bndim;j++) 
	{
	  cp[j] = bp[i][j];
	}
    by[i]=func(cp);      // Funktionswerte der Startwerte berechnen
  }
}

void
Tamoeba::ergebnis_sichern()
{
short  j;

  for (j=1;j <= bndim;j++) ap[j]=bp[1][j];
}

bool
Tamoeba::abbruch(double ftol)
{
bool  ende;
short    i;
double   fehler;

  ende = false;
  fehler = 0;
  for (i=1;i <= bndim;i++) fehler += fabs(ap[i]-bp[1][i]);
  if (fehler < ftol) ende = true;
  ergebnis_sichern();
  return ende;  
}

double 
Tamoeba::fit(double ftol, short &iter)
{
short  i;
glnp    cp;

  ende = false;
  singularMatrix = false;
  for (i=1;i < 4;i++)
  {
    startwerte();
    funktionswerte(cp);
    if (timeconstantSero) break;
    amoeba(bp,by,bndim,ftol,iter);  // neue Kanalstr�me u. Sigma berechnen
    if ((ende) || (singularMatrix) || (timeconstantSero) ||  (abbruch(ftol))) break; // !!!Reihenfolge der Abbruchbed. hier wichtig!!! warum?
  }
  if (!(Application.autom_steuer))
  {
    if (i>=3) 
    {
      warning((char *)"pause in FIT\nexecuted amoeba 3 x");
      // MessageBox(HWindow,'pause in FIT', 'executed amoeba 3 x', mb_IconExclamation);
    }
  }
  if (!(timeconstantSero))
  {
	return func(ap);
  }
  return 0; // in PASCAL vordefiniert?
}

// init und dynamisch niveaus
TGaussFit::TGaussFit(Tampl_fit_type &af)
{
short  k,i;

  n_channels = af.n_channels;
  n_extra_channels = af.n_extra_channels;
  ap[1]=af.sigma;
  ap[2]=af.i_null;
  ap[3]=af.i_channel;
  ap[4]=af.i_extra_channel;
  bndim=4;   // Dimension: sigma, i_null, i_channel, i_extra_channel
  if (ap[4] == 0) bndim = 3;
  if (ap[3] == 0) bndim = 2;

  k=determine_base();

  for (i=1;i < k+1;i++) alista[i]=i;
  ama=k;
  amfit=k;
}

void
TGaussFit::done(Tampl_fit_type &af)									
{
  af.sigma=ap[1];
  af.i_null=ap[2];
  if (bndim >= 3) af.i_channel=ap[3];
  if (bndim >= 4) af.i_extra_channel=ap[4];

  niveaus.clear();
}

// mehrere Gaussfunktionen
void 
TGaussFit::calculate(Tniveau niveau, 
                     double z, 
                     RealArrayMA afunc, 
                     short  &k, 
                     double &m, 
                     double &b)
{
  m = ap[2] + niveau.channel * ap[3] + niveau.extra_channel * ap[4];
  b = -0.5 * ((z - m) / (ap[1])) * ((z - m) / (ap[1]));
  if (fabs(b) > 50)
  {
    afunc[k] = 0;
  } else {
    afunc[k] = exp(b);
  }
  k++;
}

void 
TGaussFit::funcs(double z, RealArrayMA afunc)
{
short  k;
double m;
double b;

  // Create constant iterator for list.
  __gnu_cxx::slist<Tniveau>::const_iterator iter; //Tibias 2017

  k=1;
  // Iterate through list
  for (iter=niveaus.begin(); iter != niveaus.end(); iter++)
  {
    calculate(*iter, z, afunc, k, m, b);
  }
}

double 
TGaussFit::summe_gauss_glocken(RealArrayMA  afunc)
{
// Summe aller Gaussfunktionen

unsigned short       i;
double      ret;

  ret=0;
  for (i=1;i < niveaus.size()+1;i++)
  {
    // eigentlich d�rfte kein wert negativ sein!
    ret += ((afunc[i] * aa[i]) > 0 ? afunc[i] * aa[i] : 0);
  }
  return (ret > 0) ? ret : 0;
}        

double 
TGaussFit::gauss_glocke(short niveau, RealArrayMA  afunc)
{
// nur eine gauss glocke

  return (afunc[niveau] * aa[niveau] > 0) ? afunc[niveau] * aa[niveau] : 0;
}        

bool
TGaussFit::IsAlreadyThere(double akt)
{
// Create constant iterator for list.
__gnu_cxx::slist<Tniveau>::const_iterator iter;		//Tobias 2017

  // Iterate through list
  for (iter=niveaus.begin(); iter != niveaus.end(); iter++)
  {
    if (iter->n == akt)  return true;
  }
  return false;
}

// Bestimmung der Basisvektoren in Liste, R�ckgabe: Anz. der Vektoren
// gleiche Niveaus werden weggelassen! ->"Basis" d.h. l.u.
short  
TGaussFit::determine_base()
{ 
short  i,j;
double akt_niveau;

   niveaus.clear();
   for (i=0;i < n_channels+1;i++)
   {
     for (j=0;j < n_extra_channels+1;j++)
     {
       akt_niveau=ap[2] + (i * ap[3])+(j * ap[4]);      
       if (!(IsAlreadyThere(akt_niveau)))
       {
         Tniveau niv(akt_niveau,i,j);
         niveaus.push_front(niv);
       }
     }
   }
   return niveaus.size();
}


TExpFit::TExpFit(TTau &tau)
{
short  i;

  for (i=1;i <= NMaxTimeConstants;i++)
  {
    ap[i] = tau[i];
  }
  bndim = NMaxTimeConstants;
  for (i=NMaxTimeConstants;i >= 1;i--)
  {
    if (ap[i] == 0) bndim = bndim-1;
  }
  for (i=1;i <= bndim;i++)
  {
	alista[i] = i;
  }
  ama = bndim;   // Die Zahl der Exponentilafunktionen (ama) ist hier gleich der Zahl der fit-Parameter (bndim),
  amfit = bndim; // da jeder Exponentialfunktion eine Zeitkonstante, die zu fitten ist, enth�lt.
}

void 
TExpFit::done(TTau &tau)
{
short i;

  for (i=1;i <= bndim;i++)
  {
    tau[i] = ap[i];
  }
//  inherited done; Tamoeba hat kein done!
}

TExpFit::~TExpFit()
{
}

void 
TExpFit::funcs(double z, RealArrayMA afunc)
{
short  i;
double x;

  for (i=1;i <= bndim;i++)
  {
    if (fabs(ap[i]) < 0.000001)
    {
      timeconstantSero = true;
      ap[i] = 0.000001; // auf Wert <> Null setzen, damit kein Absturz
	}
    if (timeconstantSero)
	{
      if (!(Application.autom_steuer)) warning((char *)"Pause in FUNCS: time constant zero");
      return;
	}
    x = -z/ap[i];
    if (x > 20)
	{
      afunc[i] = exp(20);
	} else {
	  if (x < -20)
      {
        afunc[i] = exp(-20);
      } else {
        afunc[i] = exp(x);
	  }
	}
  }
}

double 
TExpFit::summe_expo(double x)
{
short  i;
double f;
RealArrayMA  afunc;

  funcs(x, afunc);
  f = 0;
  if (!timeconstantSero)
  {
    for (i=1;i <= bndim;i++) f = f + afunc[i]*aa[i];
  }
  return f;
}
