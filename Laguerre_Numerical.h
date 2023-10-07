#include<complex>
#include<vector>
#include<limits>

using namespace std;

typedef complex<double> Complex;
//typedef const vector<complex<double>> VecComplex_I;
typedef vector<complex<double>> VecComplex_I;
typedef vector<complex<double>> VecComplex, VecComplex_O;
typedef vector<complex<double>>::size_type SizeType;

void laguer(VecComplex_I &a, Complex &x, int &its){
    const int MR=8, MT=10, MAXIT=MT*MR;
    const double EPS=numeric_limits<double>::epsilon();
    // EPS here: estimated fractional roundoff error

    // we try to break (rare) limit cycles with MR 
    // MR different fractional values, once every MT steps,
    // for MAXIT total allowed iterations.
    static const double frac[MR+1]=
        {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    // Fractions used to break a limit cycle.
    Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
    int m=a.size()-1;
    for(int iter=1; iter<=MAXIT; iter++){
        // loop over iterations up to allowed maximum
        its = iter;
        b = a[m];
        double err = abs(b);
        d = f = 0.0;
        double abx = abs(x);
        for(int j=m-1; j>=0; j--){
            // efficient computation of the polynomial and
            // its first two derivatives. f stores P''/2
            f = x*f + d;
            d = x*d + b;
            b = x*b + a[j];
            err = abs(b) + abx*err;
        }
        // estimate of roundoff error in evaluating 
        // polynomial
        err *= EPS;
        if(abs(b)<=err) return; // we are on the root
        // the generic case: use Laguerre's formula
        g = d/b; 
        g2 = g*g;
        h = g2-2.0*f/b;
        sq = sqrt(((double) m-1)*(((double) m)*h-g2));
        gp = g + sq;
        gm = g - sq;
        double abp = abs(gp);
        double abm = abs(gm);
        if (abp < abm) gp = gm;
        dx = max(abp, abm) > 0.0 ? ((double) m)/gp : 
            polar(1+abx, (double) iter);
        x1 = x - dx;
        if(x == x1) return; // converged
        if(iter % MT != 0) x = x1;
        // every so often we take a fractional step, 
        // to break any limit cicle (itself a rare 
        // occurrence)
        else x -= frac[iter/MT] * dx;
    }
    throw("too many iterations in laguer");
    // very unusual: can occurr only for complex roots.
    // try a different starting guess.
}

void zroots(VecComplex_I &a, VecComplex_O &roots, const bool polish) {
    const double EPS=1.0e-14; // a small number
    int i, its;
    Complex x, b, c;
    int m = a.size() - 1;
    VecComplex ad(m+1);
    // copy of coefs. for successive deflation
    for(int j=0;j<=m;j++) ad[j] = a[j]; 
    for(int j=m-1;j>=0;j--){
        x = 0.0; // start at zero to favor convergence to
        // smallest remaining root, and return the root.
        VecComplex ad_v(j+2);
        for(int jj=0; jj<j+2; jj++) ad_v[jj] = ad[jj];
        laguer(ad_v, x, its);
        // cout << "Its: \t " << its << '\n';
        if(abs(imag(x)) <= 2.0*EPS*abs(real(x)))
            x = Complex(real(x), 0.0);
        roots[j] = x;
        b = ad[j+1];
        for(int jj=j; jj>=0; jj--){
            c = ad[jj];
            ad[jj] = b;
            b = x*b + c;
        }
    }
    if (polish)
        for(int j=1; j<m; j++){
            x = roots[j];
            for(i=j-1; j<m; j++){
                if( real(roots[i]) <= real(x)) break;
                roots[i+j] = roots[i];
            }
            roots[i+j] = x;
        }
}
