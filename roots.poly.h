#include<complex>
#include<vector>
#include<limits>
#include "ExtendedFMA.h"

using namespace std;

typedef complex<double> Complex;
typedef vector<complex<double>> VecComplex_I;
typedef vector<complex<double>> VecComplex, VecComplex_O;
typedef vector<complex<double>>::size_type SizeType;

constexpr double EPS=numeric_limits<double>::epsilon();


void laguer(VecComplex_I &a, Complex &x){
    const int MT = 10;
    static const double frac[9] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
    double err, abx, abp, abm;
    int m = a.size() - 1;
    Complex a_m = a[m];
    for(int iter=1; iter<=80; iter++){
        // loop over iterations up to allowed maximum
        b = a_m;
        err = abs(b);
        d = f = 0.0;
        abx = abs(x);
        for(int j=m-1; j>=0; j--){
            // efficient computation of the polynomial and
            // its first two derivatives. f stores P''/2
            f = implementations::fma(x, f, d);
            d = implementations::fma(x, d, b);
            b = implementations::fma(x, b, a[j]);
            err = fma(err, abx, abs(b));
        }
        // estimate of roundoff error in evaluating 
        // polynomial
        if(abs(b) <= err*EPS) return; // we are on the root
        // the generic case: use Laguerre's formula
        g = d/b; 
        g2 = pow(g, 2);
        h = implementations::fma(f/b,-static_cast<double>(2), g2);
        sq = sqrt(((double) m-1)*(implementations::fma(h,((double) m),-g2)));
        gp = g + sq;
        gm = g - sq;
        abp = abs(gp);
        abm = abs(gm);
        if (abp < abm) gp = gm;
        dx = max(abp, abm) > 0.0 ? ((double) m)/gp : polar(1+abx, (double) iter);
        x1 = x - dx;
        if(x == x1) return; // converged
        // every so often we take a fractional step, 
        // to break any limit cicle (itself a rare 
        // occurrence)
        iter % MT != 0 ? x = x1 : x = implementations::fma(x, -frac[iter/MT], dx);
    }
    throw("too many iterations in laguer");
    // very unusual: can occurr only for complex roots.
    // try a different starting guess.
}

void zroots(VecComplex_I &a, VecComplex_O &roots, const bool polish) {
    Complex x, _b, _c;
    int m = a.size() - 1;
    VecComplex ad(m+1);
    // copy of coefs. for successive deflation
    for(int j=0;j<=m;j++) ad[j] = a[j];  
    for(int j=m-1;j>=0;j--){
        x = 0.0; // start at zero to favor convergence to
        // smallest remaining root, and return the root.
        VecComplex ad_v(ad.cbegin(), ad.cbegin()+j+2); //for(int jj=0; jj<j+2; jj++) ad_v[jj] = ad[jj];
        laguer(ad_v, x);
        if(fabs(imag(x)) <= fabs(real(x))*EPS) x.imag(0);
        roots[j] = x;
        _b = ad[j+1];
        for(int jj=j; jj>=0; jj--){
            _c = ad[jj];
            ad[jj] = _b;
            _b = implementations::fma(x, _b, _c);
        }
    }
    if (polish){
        int i;
        for(int j=1; j<m; j++){
            x = roots[j];
            for(i=j-1; j<m; j++){
                if(real(roots[i]) <= real(x)) break;
                roots[i+j] = roots[i];
            }
            roots[i+j] = x;
        }
    }
}
