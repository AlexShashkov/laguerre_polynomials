#include <iostream>
#include "roots.poly.h"

int main(){
    char buffer[100];
    VecComplex_I a = {
            52 ,
            -112 ,
            24 ,
            8 ,
            1
    };
    VecComplex_O roots(a.size()-1);

    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g",real(a[i]));
        cout << "a[" << i << "]= \t" << buffer << '\n';
    }
    cout << "poly: ";
    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(a[i]),imag(a[i]));
        cout << "(" << buffer << ")" << "x^" << i;
        if(i<a.size()-1) cout << '+';
    }

    cout << '\n';

    zroots(a, roots, false);
    cout << "roots \n";
    for(SizeType i=0; i<roots.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(roots[i]),imag(roots[i]));
        cout << buffer << '\n';
    }

    cout << '\n';

    a = {
            -1000.0, 1000002.0, -2000.001, 1.0
    };
    roots = VecComplex(a.size()-1);

    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g",real(a[i]));
        cout << "a[" << i << "]= \t" << buffer << '\n';
    }
    cout << "poly: ";
    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(a[i]),imag(a[i]));
        cout << "(" << buffer << ")" << "x^" << i;
        if(i<a.size()-1) cout << '+';
    }

    cout << '\n';

    zroots(a, roots, false);
    cout << "roots \n";
    for(SizeType i=0; i<roots.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(roots[i]),imag(roots[i]));
        cout << buffer << '\n';
    }

    cout << '\n';

    a = {
            -6.855188152137764e-15,
            7.500004043464861e-01,
            2.668263562685250e+07,
            5.993293694242965e+11,
            3.621535214588474e+11,
    };
    roots = VecComplex(a.size()-1);

    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g",real(a[i]));
        cout << "a[" << i << "]= \t" << buffer << '\n';
    }
    cout << "poly: ";
    for(SizeType i=0; i<a.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(a[i]),imag(a[i]));
        cout << "(" << buffer << ")" << "x^" << i;
        if(i<a.size()-1) cout << '+';
    }

    cout << '\n';

    zroots(a, roots, false);
    cout << "roots \n";
    for(SizeType i=0; i<roots.size(); i++){
        sprintf(buffer, "%.15g+%0.15g j",real(roots[i]),imag(roots[i]));
        cout << buffer << '\n';
    }

    cout << '\n';


    return 0;
}
