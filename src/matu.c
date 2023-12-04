/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <math.h>

//#include "glib.h"

#include "matu.h"

void cdiv(double ar, double ai, double br, double bi, double &cr, double &ci)
//      complex division, (cr,ci) = (ar,ai)/(br,bi)
{
double  s,ars,ais,brs,bis;

  s = fabs(br) + fabs(bi) ;
  ars = ar / s ;
  ais = ai / s ;
  brs = br / s ;
  bis = bi / s ;
  s   = brs * brs + bis * bis ;
  cr  = (ars * brs + ais * bis) / s ;
  ci  = (ais * brs - ars * bis) / s ;
}

void 
selmhes(short n, short low, short igh, MatD a, VecID Int)
{
short  i,j,m,kp1,la,mm1;
double x,y;
//Label 100;

  la  = igh - 1;
  kp1 = low + 1;
  if (la >= kp1)
  {
    for (m=kp1;m <= la;m++)
	{
      mm1 = m - 1;
      x   = 0.0;
      i   = m ;

      for (j=m;j <= igh;j++)
	  {
        if (fabs(a[j][mm1]) > fabs(x))
		{
          x = a[j][mm1];
          i = j;
		}
	  }
      Int[m] = i;

      if (i != m)
	  {
        //     .......... interchange rows and columns of a ..........
        for (j=mm1;j <= n;j++) 
		{
          y = a[i][j];
          a[i][j] = a[m][j];
          a[m][j] = y;
		}
        for (j=1;j <= igh;j++) 
		{
          y = a[j][i];
          a[j][i] = a[j][m];
          a[j][m] = y;
		}
        //     .......... end interchange ..........
	  }
      if (x != 0.0)
	  {
        for (i=m+1;i <= igh;i++)
		{
          y = a[i][mm1];
          if (y != 0.0)
		  {
            y = y / x ;
            a[i][mm1] = y ;
            for (j=m;j <= n;j++) a[i][j] -= y * a[m][j];
            for (j=1;j <= igh;j++) a[j][m] += y * a[j][i];
          }
        }
      }
    }
  }
}

void
seltran(short n, short low, short igh, MatD a, VecID Int, MatD z)
{
short  i,j,mp,mp1;
//    .......... initialize z to identity matrix ..........
//  fillchar (z,sizeof(z),#0) ;
// F�llt einen Speicherbereich von der Gr��e Count mit dem Wert value (vom Typ Byte oder Char).

  for (i=1;i <= n;i++) z[i][i] = 1;

  if (igh - low >= 2) 
  {
    for (mp=igh-1;mp >= low+1;mp--)
	{
      mp1 = mp + 1;
      for (i=mp1;i <= igh;i++) z[i][mp] = a[i][mp-1];
      i = Int[mp];
      if (i != mp)
	  {
       for (j=mp;j <= igh;j++)
	   {
         z[mp][j] = z[i][j];
         z[i][j] = 0;
       }
       z[i][mp] = 1;
	  }
    }
  }
}

void 
selmbak(short low, short igh, short m, MatD a, VecID Int, MatD z)
{
short  i,j,la,kp1,mp1;
unsigned short mp;
double x;

  la = igh - 1;
  kp1 = low + 1;
  if ((m > 0) && (la >= kp1))
  {
    for (mp=la;mp >= kp1;mp--)
	{
      mp1 = mp + 1;
      for (i=mp1;i <= igh;i++)
	  {
        x = a[i][mp-1];
        if (x != 0)
		{
          for (j=1;j <= m;j++) z[i][j] += x * z[mp][j];
		}
      }
      i = Int[mp];
      if (i != mp)
	  {
        for (j=1;j <= m;j++)
		{
          x = z[i][j];
          z[i][j] = z[mp][j];
          z[mp][j] = x;
        }
      }
    }
  }
}

void 
sbalanc(short n, short &low, short &igh, MatD a, VecD scale)
{
short   i,j,k,l;
double  c,f,g,r,s,b2,radix;
bool noconv;

  radix = 2.0e0; // 16.0e0 ; fuer IBM Mainframe
  b2 = radix * radix;
  k = 1;
  l = n;
  j = l+1;
  do
  {
    j--;
    i = 0;
    do
	{
      i++;
      if ((i != j) && (a[j][i] != 0)) i = l + 1;
    } while (i < l) ;
    if (i ==l)
	{
      scale[l] = j;
      if (j != l)
	  {
        //     .......... procedure for row and column exchange .......... 
        for (i=1;i <= l;i++)
		{
         f = a[i][j];
         a[i][j] = a[i][l];
         a[i][l] = f;
        }
        for (i=k;i <= n;i++)
		{
          f = a[j][i];
          a[j][i] = a[l][i];
          a[l][i] = f;
		}        
	  }
      if (l == 1)
	  {
        low = 1;
        igh = 1;
        return;
      }
      l--;
      j = l;
    }
  } while (j != 1);
  //     .......... search for columns isolating an eigenvalue
  //                 and push them left ..........
  j = k - 1;
  do
  {
    j++;
    i = k - 1;
    do
	{
      i++ ;
      if ((i != j) && (a[j][i] != 0)) i = l+1;
    } while (i < l);
    if (i == l)
	{
      scale[k] = j;
      if (j != k)
	  {
        //     .......... procedure for row and column exchange .......... 
        for (i=1;i <= l;i++)
		{
          f = a[i][j];
          a[i][j] = a[i][k];
          a[i][k] = f;
        }
        for (i=k;i <= n;i++)
		{
          f = a[j][i];
          a[j][i] = a[k][i];
          a[k][i] = f;
        }
	  }
      k++; //naechstes K
      j = k;
    }
  } while (j != l);
  //     .......... now balance the submatrix in rows k to l ..........
  for (i=k;i <= l;i++) scale[i] = 1;
  //     .......... iterative loop for norm reduction ..........

  do
  {
    noconv = false;
    for (i=k;i <= l;i++)
	{
      c = 0;
      r = c;
      for (j=k;j <= l;j++)
	  {
        if (j != i)
		{
          c += fabs(a[j][i]);
          r += fabs(a[i][j]);
		}
	  }
      //     .......... guard against zero c or r due to underflow ..........
      if ((c != 0) && (r != 0))
	  {
        g = r / radix;
        f = 1;
        s = c + r;
        while (c < g)
		{
          f = f * radix;
          c = c * b2;
		}
        g = r * radix;
        while (c >= g)
		{
          f = f / radix;
          c = c / b2;
		}
        //     .......... now balance ..........
        if ((c + r) / f < 0.95e0 * s)
		{
          g = 1.0e0 / f;
          scale[i] = scale[i] * f;
          noconv = true;
          for (j=k;j <= n;j++) a[i][j] *= g;
          for (j=1;j <= l;j++) a[j][i] *= f;
		}
	  }
	}
  } while (noconv);
  low = k;
  igh = l;
}


void 
sbalbak(short n, short low, short igh, short m, VecID scale, MatD z)
{
short  i,j,k,ii;
double s;

  if (m > 0) 
  {
    if (igh > low)
	{
      for (i=low;i <= igh;i++)
	  {
        s = scale[i] ;
        // .... left hand eigenvectors are back transformed
        //      if the foregoing statement is replaced by
        //      s:=1.0e0/scale(i). ..........
        for (j=1;j <= m;j++) z[i][j] *= s;
      }
      for (ii=1;ii <= n;ii++)
	  {
        i = ii;
        if ((i < low) || (i > igh))
		{
          if (i < low) i = low - ii;
          k = scale[i];
          if (k != i)
		  {
            for (j=1;j <= m;j++)
			{
              s = z[i][j];
              z[i][j] = z[k][j];
              z[k][j] = s;
			}
          }
		}
	  }
    }
  }
}

void 
shqr2(short n, short low, short igh, MatD h, VECn wr, VECn wi, MatD z, short &ierr)
{
short   i,j,k,l,m,en,na,itn,its,mp2,enm2,root;
double  p=0,q=0,r=0,s=0,t=0,w=0,x=0,y=0,ra=0,sa=0,vi=0,vr=0,zz=0,norm=0,tst1=0,tst2=0;
bool notlas;

  ierr = 0;
  //  .......... store roots isolated by balanc
  //             and compute matrix norm ..........
  norm = 0;
  k = 1;
  for (i=1;i <= n;i++)
  {
    for (j=k;j <= n;j++)
	{
      norm += fabs(h[i][j]);
    }
    k = i;
    if ((i < low) || (i > igh))
	{
      wr[i] = h[i][i];
      wi[i] = 0;
    }
  }

  en = igh;
  t = 0;
  itn = 30 * n;

  //     .......... search for next eigenvalues ..........
  while (en >= low)
  {
    its = 0;
    na = en - 1;
    enm2 = na - 1;

  //     .......... look for single small sub-diagonal element
  //                for l=en step -1 until low do -- ..........
    do
	{
      l = en + 1;
      do
	  {
        l--;
        if (l > low)
		{
          s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
          if (s == 0) s = norm;
          tst1 = s;
          tst2 = tst1 + fabs(h[l][l-1]);
        }
      } while (!((tst2 == tst1) || (l == low)));

      // Der n�chste Abschnitt bis 100 wurde eingeklammert, ich wei� nicht warum!!! 
       
      // for l := en downto low do begin
      //  if (l>low) then begin
      //   s := abs(h[l-1,l-1]) + abs(h[l,l]) ;
      //   if (s = 0.0e0) then s := norm ;
      //   tst1 := s ;
      //   tst2 := tst1 + abs(h[l,l-1]) ;
      //   if (tst2 = tst1) then Goto 100 ;
      //  end ;
      // end ;
      //100 :

      //     .......... form shift .......... 
      root = 0;
      x = h[en][en];
      if (l == en)
      {
        root = 1; //one root found , 1*
      } else {
        y = h[na][na];
        w = h[en][na] * h[na][en];
      }
      if (l == na) root = 2; //two roots found , 1*

      switch (root)
	  {
        case 0:
          if (itn == 0)
		  {
            //  ... set error -- all eigenvalues have not
            //      converged after 30*n iterations ..........
            ierr = en;
            return;
		  }

          if ((its == 10) || (its == 20))
		  {
            //     .......... form exceptional shift .......... 
            t += x;
            for (i=low;i <= en;i++) h[i][i] -= x;
            s = fabs(h[en][na]) + fabs(h[na][enm2]);
            x = 0.75e0 * s;
            y = x;
            w = -0.4375e0 * s * s ;
          }
          its = its + 1;
          itn = itn - 1;

          //  ... look for two consecutive small
          //      sub-diagonal elements.
          m = enm2;
          for (m=enm2;m >= l;m--)
		  {
            zz = h[m][m];
            r = x - zz;
            s = y - zz;
            p = (r * s - w) / h[m+1][m] + h[m][m+1]; // problem klammern ?
            q = h[m+1][m+1] - zz - r - s;
            r = h[m+2][m+1];
            s = fabs(p) + fabs(q) + fabs(r);
            p = p / s;
            q = q / s;
            r = r / s;
            if (m != l) 
			{
              tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
              tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
              if (tst2 == tst1) goto label150;
            }
          }
          label150:
          m=m+1; //C geht bei der For Schleife einen zu weit runter, anders als bei Pascal
                 // hier der Ausgleich

          mp2 = m + 2;
  
          for (i=mp2;i <= en;i++)
		  {
            h[i][i-2] = 0;
            if (i != mp2) h[i][i-3] = 0;
          }

          //  ... double qr step involving rows l to en and
          //      columns m to en ..........
          for (k=m;k <= na;k++)
		  {
            notlas = (k != na);
            if (k != m)
			{
              p = h[k][k-1];
              q = h[k+1][k-1];
              r = 0;
              if (notlas) r = h[k+2][k-1];
              x = fabs(p) + fabs(q) + fabs(r);
              if (x == 0) goto label260;
              p = p / x ;
              q = q / x ;
              r = r / x ;
            }
            s = sqrt(p * p + q * q + r * r);
            if (p < 0) s = -s;

            if (k != m) 
			{
			  h[k][k-1] = -s * x;
			} else {
			  if (l != m) h[k][k-1] = -h[k][k-1];
			}
            p  = p + s;
            x  = p / s;
            y  = q / s;
            zz = r / s;
            q  = q / p;
            r  = r / p;
            if (!notlas)
			{
              //   ... row modification .......... 
              for (j=k;j <= n;j++)
			  {
                p = h[k][j] + q * h[k+1][j];
                h[k][j] = h[k][j] - p * x;
                h[k+1][j] = h[k+1][j] - p * y;
			  }
              j = en;
              if (k+3 < j) j = k+3;
              //   ... column modification .......... 
              for (i=1;i <= j;i++)
			  {
                p = x * h[i][k] + y * h[i][k+1];
                h[i][k] = h[i][k] - p;
                h[i][k+1] = h[i][k+1] - p * q;
              }
              //   ... accumulate transformations ..........
              for (i=low;i <= igh;i++)
			  {
                p = x * z[i][k] + y * z[i][k+1];
                z[i][k] = z[i][k] - p;
                z[i][k+1] = z[i][k+1] - p * q;
              }
            } else { //notlas = true
              //   ... row modification .......... 
              for (j=k;j <= n;j++)
			  {
                p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
                h[k][j] = h[k][j] - p * x;
                h[k+1][j] = h[k+1][j] - p * y;
                h[k+2][j] = h[k+2][j] - p * zz;
              }
              j = en;
              if (k+3 < j) j = k+3;

              //   ... column modification .......... 
              for (i=1;i <= j;i++)
			  {
                p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
                h[i][k] = h[i][k] - p;
                h[i][k+1] = h[i][k+1] - p * q;
                h[i][k+2] = h[i][k+2] - p * r;
              }
              //   ... accumulate transformations .......... 
              for (i=low;i <= igh;i++)
			  {
                p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
                z[i][k] = z[i][k] - p;
                z[i][k+1] = z[i][k+1] - p * q;
                z[i][k+2] = z[i][k+2] - p * r;
              }
            }
            label260: ;
          }
          break;
        
        case 1:
          //   ... one root found .......... 
          h[en][en] = x + t;
          wr[en] = h[en][en];
          wi[en] = 0;
          en = na; // eine Weniger
          break;

        case 2:
          //   ... two roots found .......... 
          p = (y - x) / 2.0e0;
          q = p * p + w;
          zz = sqrt(fabs(q));
          h[en][en] = x + t;
          x = h[en][en];
          h[na][na] = y + t;
          if (q >= 0)
		  {
            //   ... real pair .......... 
            if (p < 0) zz = -zz;
            zz = p + zz;
            wr[na] = x + zz;
            wr[en] = wr[na];
            if (zz != 0.0e0) wr[en] = x - w / zz;
            wi[na] = 0.0e0;
            wi[en] = 0.0e0;
            x = h[en][na];
            s = fabs(x) + fabs(zz);
            p = x / s;
            q = zz / s;
            r = sqrt(p * p + q * q);
            p = p / r;
            q = q / r;
            //   ... row modification ..........
            for (j=na;j <= n;j++)
			{
              zz = h[na][j];
              h[na][j] = q * zz + p * h[en][j];
              h[en][j] = q * h[en][j] - p * zz;
            }
            //   ... column modification ..........
            for (i=1;i <= en;i++)
			{
              zz = h[i][na];
              h[i][na] = q * zz + p * h[i][en];
              h[i][en] = q * h[i][en] - p * zz;
            }
            //   ... accumulate transformations .......... 
            for (i=low;i <= igh;i++)
			{
              zz = z[i][na];
              z[i][na] = q * zz + p * z[i][en];
              z[i][en] = q * z[i][en] - p * zz;
            }
          } else {
            //     .......... complex pair .......... 
            wr[na] = x + p;
            wr[en] = x + p;
            wi[na] = zz;
            wi[en] = -zz;
          }
          en = enm2; // zwei weniger
          break; // of root = 2
		default:
		  // kann nicht erreicht werden
		  break;
      }
    } while (root <= 0);
  }
  // Teil II :
  // .... all roots found.  backsubstitute to find
  //      vectors of upper triangular form .......

  if (norm == 0) return;
  for (en=n;en >= 1;en--)
  {
    p = wr[en];
    q = wi[en];
    na = en - 1;
    if (q == 0)
	{
	  // ... real vector ...
      m = en;
      h[en][en] = 1.0e0;
      if (na != 0)
	  {
        for (i=na;i >= 1;i--)
		{
          w = h[i][i] - p;
          r = 0;
          for (j=m;j <= en;j++) r += h[i][j] * h[j][en];

          if (wi[i] < 0)
		  {
            zz = w;
            s = r;
          } else {
            m = i;
            if (wi[i] == 0)
			{
              t = w;
              if (t == 0)
			  {
                tst1 = norm;
                t = tst1;
                do
				{
                  t = 0.01e0 * t;
                  tst2 = norm + t;
                } while (tst2 > tst1);
              }
              h[i][en] = -r / t;
            } else { // wi[i] > 0
			  //     .......... solve real equations ..
              x = h[i][i+1];
              y = h[i+1][i];
              q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
              t = (x * s - zz * r) / q;
              h[i][en] = t;
              if (fabs(x) > fabs(zz)) 
			  {
				h[i+1][en] = (-r - w * t) / x;
			  } else {
				h[i+1][en] = (-s - y * t) / zz;
			  }
            }
			// ... overflow control ...
            t = fabs(h[i][en]);
            if (t != 0)
			{
              tst1 = t;
              tst2 = tst1 + 1.0e0 / tst1;
              if (tst2 <= tst1)
			  {
                for (j=i;j <= en;j++) h[j][en] = h[j][en] / t;
              }
            }
          }
        }
	  }
    } else {     // ... end real vector ... 
                 // Blockende
      if (q < 0)
	  {
	    // ... begin complex vector ...
        m = na;
		// ... last vector component chosen imaginary so that
        //     eigenvector matrix is triangular ...
        if (fabs(h[en][na]) > fabs(h[na][en]))
		{
          h[na][na] = q / h[en][na];
          h[na][en] = -(h[en][en] - p) / h[en][na];
        } else {
          cdiv(0, -h[na][en], h[na][na]-p, q, h[na][na], h[na][en]);
        }
        h[en][na] = 0.0e0;
        h[en][en] = 1.0e0;
        enm2 = na - 1;
        if (enm2 != 0) 
		{
          for (i=enm2;i >= 1;i--)
		  {
            w = h[i][i] - p;
            ra = 0.0e0;
            sa = 0.0e0;

            for (j=m;j <= en;j++)
			{
              ra += h[i][j] * h[j][na];
              sa += h[i][j] * h[j][en];
            }

            if (wi[i] < 0.0e0)
			{
              zz = w;
              r = ra;
              s = sa;
            } else {
              m = i;
              if (wi[i] == 0.0e0)
			  {
                cdiv(-ra, -sa, w, q, h[i][na], h[i][en]);
              } else {
			    // ... solve complex equations ...
                x = h[i][i+1];
                y = h[i+1][i];
                vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
                vi = (wr[i] - p) * 2.0e0 * q;
                if ((vr == 0.0e0) && (vi == 0.0e0))
				{
                  tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
                  vr = tst1;
                  do
				  {
                    vr = 0.01e0 * vr;
                    tst2 = tst1 + vr;
                  } while (tst2 > tst1);
                }
                cdiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, h[i][na], h[i][en]);
                if (fabs(x) > fabs(zz) + fabs(q))
				{
                  h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
                  h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
                } else {
                  cdiv(-r - y * h[i][na], -s - y * h[i][en], zz, q, h[i+1][na], h[i+1][en]);
                }
              } 
			  // ... overflow control ...
              t = fabs(h[i][na]);
              if (t < fabs(h[i][en])) t = fabs(h[i][en]);
              if (t != 0.0e0)
			  {
                tst1 = t;
                tst2 = tst1 + 1.0e0 / tst1;
                if (tst2 <= tst1)
				{
                  for (j=i;j <= en;j++)
				  {
                    h[j][na] = h[j][na] / t;
                    h[j][en] = h[j][en] / t;
                  }
                }
              }
            } //of else
          } //for i
		} //enm2 <> 0
	  } //of q<0.0
	}
	//     .......... end complex vector ..........
  }
  //     .......... end back substitution.
  //           vectors of isolated roots ..........
  for (i=1;i <= n;i++)
  {
    if ((i < low) || (i > igh))
	{
       for (j=i;j <= n;j++) z[i][j] = h[i][j];
	}
  }
  // ... multiply by transformation matrix to give
  //     vectors of original full matrix.
  for (j=n;j >= low;j--)
  {
    m = j;
    if (igh < m) m = igh;
    for (i=low;i <= igh;i++)
	{
      zz = 0.0e0;
      for (k=low;k <= m;k++) zz += z[i][k] * h[k][j];
      z[i][j] = zz;
    }
  }
}
