#ifndef _MFA_EX_FNS
#define _MFA_EX_FNS

#include <set>
#include <mfa/types.hpp>
#include "domain_args.hpp"

// Define list of example keywords
set<string> analytical_signals = {"sine", "cosine", "sinc", "psinc1", "psinc2", "psinc3", "psinc4", "ml", "f16", "f17", "f18","sinc_sum","sinc_sum_2","rastrigin","gaussian_pair1","gaussian_pair2","ackley","schwefel","gaussian_mixture","rotating_gaussian"};
set<string> datasets_4d = {"tornado4d"};
set<string> datasets_3d = {"s3d", "nek", "rti", "miranda", "tornado"};
set<string> datasets_2d = {"cesm","vortex_street","boussinesq","hurricane_isabel","s3d"};
set<string> datasets_unstructured = {"edelta", "climate", "nuclear", "nasa"};

// REMOVE:
// Computes the analytical line integral from p1 to p2 of sin(x)sin(y). Use for testing
template<typename T>
T sintest( const VectorX<T>& p1,
           const VectorX<T>& p2)
{
    T x1 = p1(0);
    T y1 = p1(1);
    T x2 = p2(0);
    T y2 = p2(1);

    T mx = x2 - x1;
    T my = y2 - y1;
    T fctr = sqrt(mx*mx + my*my);

    // a1 = mx, b1 = x1, a2 = my, b2 = y1
    T int0 = 0, int1 = -7;
    if (mx != my)
    {
        int0 = 0.5*(sin(x1-y1+(mx-my)*0)/(mx-my) - sin(x1+y1+(mx+my)*0)/(mx+my));
        int1 = 0.5*(sin(x1-y1+(mx-my)*1)/(mx-my) - sin(x1+y1+(mx+my)*1)/(mx+my));
    }
    else
    {
        int0 = 0.5 * (0*cos(x1-y1) + sin(x1+y1+(mx+my)*0)/(mx+my));
        int1 = 0.5 * (1*cos(x1-y1) + sin(x1+y1+(mx+my)*1)/(mx+my));
    }
    
    return (int1 - int0) * fctr;
}

// evaluate sine function
template<typename T>
void sine(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    T retval = 1.0;
    for (auto i = 0; i < domain_pt.size(); i++)
        retval *= sin(domain_pt(i) * args.f[k]);
    retval *= args.s[k];

    for (int l = 0; l < output_pt.size(); l++)
    {
        output_pt(l) = retval * (1+l);
    }

    return;
}

template<typename T>
void cosine(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    T retval = 1.0;
    for (auto i = 0; i < domain_pt.size(); i++)
        retval *= cos(domain_pt(i) * args.f[k]);
    retval *= args.s[k];

    for (int l = 0; l < output_pt.size(); l++)
    {
        output_pt(l) = retval * (1+l);
    }

    return;        
}

// evaluate the "negative cosine plus one" (f(x) = -cos(x)+1) function
// used primarily to test integration of sine
template<typename T>
void ncosp1(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs& args, int k)
{
    if (output_pt.size() != 1)
    {
        fprintf(stderr, "Error: ncosp1 only defined for scalar output.\n");
        exit(0);
    }

    T retval = 1.0;
    for (auto i = 0; i < domain_pt.size(); i++)
        retval *= 1 - cos(domain_pt(i) * args.f[k]);
    retval *= args.s[k];

    output_pt(0) = retval;
    return;        
}

// evaluate sinc function
template<typename T>
void sinc(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs& args, int k)
{
    T retval = 1.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        if (domain_pt(i) != 0.0)
            retval *= (sin(domain_pt(i) * args.f[k] ) / domain_pt(i));
    }
    retval *= args.s[k];

    for (int l = 0; l < output_pt.size(); l++)
    {
        output_pt(l) = retval * (1+l);
    }

    return;
}

//sin(5x)/(5x)+sin(5y)/(5y)+sin(5z)/(5z)
template<typename T>
void sinc_sum(const VectorX<T>&  domain_pt,VectorX<T>& output_pt,
        DomainArgs&  args)
{
    T coe = 5;
    DomainArgs* a = &args;

    T retval = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        if (domain_pt(i) != 0.0)
        {

            retval += 5.0*sin(coe * domain_pt(i)) / (coe *domain_pt(i));
            
        }
        else
        {
            retval += 5.0;
        }               
    }
    retval *= a->s[0];

    output_pt(0) = retval;
    return;
}


    //sin(3x+2)/(3x+2)+sin(3y+2)/(3y+2)+sin(3z+2)/(3z+2)
template<typename T>
void sinc_sum_2(const VectorX<T>&  domain_pt,VectorX<T>& output_pt,
        DomainArgs&  args)
{
    T coe = 3;
    DomainArgs* a = &args;

    T retval = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        if ((coe *domain_pt(i)+2.0) != 0.0)
        {

            retval += 5.0*sin(coe * domain_pt(i)+2.0) / (coe *domain_pt(i)+2.0);
            
        }
        else
        {
            retval += 5.0;
        }               
    }
    retval *= a->s[0];

    output_pt(0) = retval;

    return;
}


//ackley function https://www.sfu.ca/~ssurjano/ackley.html
template<typename T>
void ackley(const VectorX<T>&   domain_pt,VectorX<T>& output_pt,DomainArgs& args)
{
    T a=20;
    T b=0.2;
    T c=2.0*M_PI;
    T d=domain_pt.size();

    T term_0=-b*sqrt(domain_pt.squaredNorm()/d);
    T term_1=0.0;
    for(int i=0;i<domain_pt.size();i++)
    {
        term_1+=cos(c*domain_pt[i]);
    }
    
    output_pt(0)= args.s[0] *(-a*exp(term_0)-exp(term_1/d)+a+M_E);
    return;

}



// only for 2d domain
template<typename T>
void gaussian_mixture(const VectorX<T>&   domain_pt,VectorX<T>& output_pt)
{
    T x = domain_pt(0);
    T y = domain_pt(1);

    T f = exp(-(8.0*(x + 0.4)*(x + 0.4) + 4*y*y)) + exp(-8.0*(x-0.5)*(x-0.5)-4.0*y*y)+exp(-8.0*x*x-4.0*(y-0.77)*(y-0.77))+exp(-8.0*x*x-4.0*(y-1.5)*(y-1.5))+0.2*exp(-0.3*x*x-0.3*(y-0.5)*(y-0.5));
    
    output_pt(0) = f;

    return;

}


// template<typename T>
// void rotating_gaussian(const VectorX<T>&   domain_pt,VectorX<T>& output_pt)
// {
//     T cx = 0.0;
//     T cy = 0.0;
//     T r1 = 0.7; T r2 = 0.7;
//     T phi1 = 0; T phi2 = M_PI; 
//     T omega1 = 0.5 * M_PI; T omega2 = 0.5 * M_PI;
//     T x1= cx + r1 * cos(omega1 * domain_pt(2) + phi1);
//     T y1 = cy + r1 * sin(omega1 * domain_pt(2) + phi1);
    
//     T x2 = cx + r2 * cos(omega2 * domain_pt(2) + phi2);
//     T y2 = cy + r2 * sin(omega2 * domain_pt(2) + phi2);

//     T sigmax1 = 0.5; T sigmay1 = 0.5;
//     T sigmax2 = 0.4; T sigmay2 = 0.4;
//     T A1 = 1.0; T A2 = 0.8;

//     T f1 = A1 * exp(-((domain_pt(0) - x1)*(domain_pt(0) - x1)/(2*sigmax1*sigmax1) + (domain_pt(1) - y1)*(domain_pt(1) - y1)/(2*sigmay1*sigmay1)));
//     T f2 = A2 * exp(-((domain_pt(0) - x2)*(domain_pt(0) - x2)/(2*sigmax2*sigmax2) + (domain_pt(1) - y2)*(domain_pt(1) - y2)/(2*sigmay2*sigmay2)));
//     output_pt(0) = f1 + f2;

//     return;
// }
template<typename T>
void rotating_gaussian(const VectorX<T>&   domain_pt,VectorX<T>& output_pt)
{
    T cx = 0.0;
    T cy = 0.0;
    T r1 = 0.7; 
    T phi1 = 0; 
    T omega1 = 0.5 * M_PI;
    T sigmax1 = 0.5; T sigmay1 = 0.5;
    T sigmax2 = 0.5; T sigmay2 = 0.5;
    T A1 = 1.0; T A2 = 1.0;

    T x1= cx + r1 * cos(omega1 * domain_pt(2) + phi1);
    T y1 = cy + r1 * sin(omega1 * domain_pt(2) + phi1);

    T f2=0.0;

    T f1 = A1 * exp(-((domain_pt(0) - x1)*(domain_pt(0) - x1)/(2*sigmax1*sigmax1) + (domain_pt(1) - y1)*(domain_pt(1) - y1)/(2*sigmay1*sigmay1)));


    if (domain_pt(2)<3.0 && domain_pt(2)>1.0)
    {
        // x1*=1.25-0.25*domain_pt(2);
        // y1*=1.25-0.25*domain_pt(2);

        T x2= cx + r1 * cos(-omega1 * domain_pt(2) + M_PI + phi1);
        T y2 = cy + r1 * sin(-omega1 * domain_pt(2) + M_PI + phi1);

        f2 = A2 * exp(-((domain_pt(0) - x2)*(domain_pt(0) - x2)/(2*sigmax2*sigmax2) + (domain_pt(1) - y2)*(domain_pt(1) - y2)/(2*sigmay2*sigmay2)));

        if(f1 < f2)
        {
            f1 = f2;
        }

    }

    output_pt(0) = f1;

    return;
}

template<typename T>
void gaussian_pair1(const VectorX<T>&   domain_pt,VectorX<T>& output_pt)
{
    T x = domain_pt(0);
    T y = domain_pt(1);

    T x0=0.5;
    T y0=0.4;
    T sigma=0.1;
    T f = 0.25*exp(-(x-x0)*(x-x0)/(2*sigma*sigma)-(y-y0)*(y-y0)/(2*sigma*sigma));
    output_pt(0) = f;
    return;

}


template<typename T>
void gaussian_pair2(const VectorX<T>&   domain_pt,VectorX<T>& output_pt)
{
    T x = domain_pt(0);
    T y = domain_pt(1);

    T x0=0.3;
    T y0=0.2;
    T sigma=0.1;

    T x1=0.75;
    T y1=0.25;
    T sigma2=0.12;

    T f = 0.25*exp(-(x-x0)*(x-x0)/(2*sigma*sigma)-(y-y0)*(y-y0)/(2*sigma*sigma))
    +0.5*exp(-(x-x1)*(x-x1)/(2*sigma2*sigma2)-(y-y1)*(y-y1)/(2*sigma2*sigma2));

    output_pt(0) = f;
    return;

}


template<typename T>
void schwefel(const VectorX<T>&domain_pt,VectorX<T>& output_pt, DomainArgs& args)
{
    T a=418.9829;
    T d=domain_pt.size();

    T term_1=0.0;
    for(int i=0;i<domain_pt.size();i++)
    {
        term_1+=domain_pt[i]*sin(sqrt(abs(domain_pt[i])));
    }
    
    output_pt(0) = args.s[0]*0.5 *(a*d-term_1);
    return;

}

template<typename T>
void rastrigin(const VectorX<T>&   domain_pt,VectorX<T>& output_pt, DomainArgs& args)
{
    T a=10;
    T retval = a*domain_pt.size();
    for(int i=0;i<domain_pt.size();i++)
    {
        retval += domain_pt[i]*domain_pt[i]-a*cos(2.0*M_PI*domain_pt[i]);
    }

    retval *= 0.1;

    output_pt(0)=args.s[0] *retval;
    return;
}

// evaluate n-d poly-sinc function version 1
template<typename T>
void polysinc1(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    int dim = output_pt.size();

    // a = (x + 1)^2 + (y - 1)^2 + (z + 1)^2 + ...
    // b = (x - 1)^2 + (y + 1)^2 + (z - 1)^2 + ...
    T a = 0.0;
    T b = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        T s, r;
        if (i % 2 == 0)
        {
            s = domain_pt(i) + 1.0;
            r = domain_pt(i) - 1.0;
        }
        else
        {
            s = domain_pt(i) - 1.0;
            r = domain_pt(i) + 1.0;
        }
        a += (s * s);
        b += (r * r);
    }

    // Modulate sinc shape for each output coordinate if science variable is vector-valued
    for (int l = 0; l < dim; l++)
    {
        // a1 = sinc(a*(1+l)); b1 = sinc(b*(1+l))
        T a1 = (a == 0.0 ? 1.0*(1+l) : sin(a*(1+l)) / a);
        T b1 = (b == 0.0 ? 1.0*(1+l) : sin(b*(1+l)) / b);
        output_pt(l) = args.s[k] * (a1 + b1);               // scale entire science variable by s[k]
    }

    return;
}

// evaluate n-d poly-sinc function version 2
template<typename T>
void polysinc2(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    int dim = output_pt.size();

    // a = x^2 + y^2 + z^2 + ...
    // b = 2(x - 2)^2 + (y + 2)^2 + (z - 2)^2 + ...
    T a = 0.0;
    T b = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        T s, r;
        s = domain_pt(i);
        if (i % 2 == 0)
        {
            r = domain_pt(i) - 2.0;
        }
        else
        {
            r = domain_pt(i) + 2.0;
        }
        a += (s * s);
        if (i == 0)
            b += (2.0 * r * r);
        else
            b += (r * r);
    }

    // Modulate sinc shape for each output coordinate if science variable is vector-valued
    for (int l = 0; l < dim; l++)
    {
        // a1 = sinc(a*(1+l)); b1 = sinc(b*(1+l))
        T a1 = (a == 0.0 ? 1.0*(1+l) : sin(a*(1+l)) / a);
        T b1 = (b == 0.0 ? 1.0*(1+l) : sin(b*(1+l)) / b);
        output_pt(l) = args.s[k] * (a1 + b1);               // scale entire science variable by s[k]
    }

    return;
}

// evaluate n-d poly-sinc function version 3 (differs from version 2 in definition of a)
template<typename T>
void polysinc3(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    int dim = output_pt.size();

    // a = sqrt(x^2 + y^2 + z^2 + ...)
    // b = 2(x - 2)^2 + (y + 2)^2 + (z - 2)^2 + ...
    T a = 0.0;
    T b = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        T s, r;
        s = domain_pt(i);
        if (i % 2 == 0)
        {
            r = domain_pt(i) - 2.0;
        }
        else
        {
            r = domain_pt(i) + 2.0;
        }
        a += (s * s);
        if (i == 0)
            b += (2.0 * r * r);
        else
            b += (r * r);
    }
    a = sqrt(a);

    // Modulate sinc shape for each output coordinate if science variable is vector-valued
    for (int l = 0; l < dim; l++)
    {
        // a1 = sinc(a*(1+l)); b1 = sinc(b*(1+l))
        T a1 = (a == 0.0 ? 1.0*(1+l) : sin(a*(1+l)) / a);
        T b1 = (b == 0.0 ? 1.0*(1+l) : sin(b*(1+l)) / b);
        output_pt(l) = args.s[k] * (a1 + b1);               // scale entire science variable by s[k]
    }

    return;
}

// evaluate n-d poly-sinc function version 4 (differs from version 3 in suppresion of signal at larger radii)
template<typename T>
void polysinc4(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    int dim = output_pt.size();

    // a = sqrt(x^2 + y^2 + z^2 + ...)
    // b = 2(x - 2)^2 + (y + 2)^2 + (z - 2)^2 + ...
    // c = (0.1 + x^2 + y^2 + z^2 + ...)^1/2
    T a = 0.0;
    T b = 0.0;
    T c = 0.1;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        T s, r;
        s = domain_pt(i);
        if (i % 2 == 0)
        {
            r = domain_pt(i) - 2.0;
        }
        else
        {
            r = domain_pt(i) + 2.0;
        }
        a += (s * s);
        if (i == 0)
            b += (2.0 * r * r);
        else
            b += (r * r);
        c += domain_pt(i) * domain_pt(i);
    }
    a = sqrt(a);
    c = pow(c, 0.5);

    // Modulate sinc shape for each output coordinate if science variable is vector-valued
    for (int l = 0; l < dim; l++)
    {
        // a1 = sinc(a); b1 = sinc(b)
        T a1 = (a == 0.0 ? 1.0 : sin(a) / a);
        T b1 = (b == 0.0 ? 1.0 : sin(b) / b);
        output_pt(l) = args.s[k] * (a1 + b1) / c;       // scale entire science variable by s[k]
    }

    return;
}

// evaluate Marschner-Lobb function [Marschner and Lobb, IEEE VIS, 1994]
// only for a 3d domain
// using args f[0] and s[0] for f_M and alpha, respectively, in the paper
template<typename T>
void ml(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    if (domain_pt.size() != 3 || output_pt.size() != 1)
    {
        fprintf(stderr, "Error: Marschner-Lobb function is only defined for a 3d domain and scalar output.\n");
        exit(0);
    }
    T fm           = args.f[0];
    T alpha        = args.s[0];
//         T fm = 6.0;
//         T alpha = 0.25;
    T x            = domain_pt(0);
    T y            = domain_pt(1);
    T z            = domain_pt(2);

    T rad       = sqrt(x * x + y * y + z * z);
    T rho       = cos(2 * M_PI * fm * cos(M_PI * rad / 2.0));
    T retval    = (1.0 - sin(M_PI * z / 2.0) + alpha * (1.0 + rho * sqrt(x * x + y * y))) / (2 * (1.0 + alpha));

    output_pt(0) = retval;
    return;
}

// evaluate f16 function
template<typename T>
void f16(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    if (domain_pt.size() != 2 || output_pt.size() != 1)
    {
        fprintf(stderr, "Error: f16 function is only defined for a 2d domain and scalar output.\n");
        exit(0);
    }

    T retval =
        (pow(domain_pt(0), 4)                        +
        pow(domain_pt(1), 4)                        +
        pow(domain_pt(0), 2) * pow(domain_pt(1), 2) +
        domain_pt(0) * domain_pt(1)                 ) /
        (pow(domain_pt(0), 3)                        +
        pow(domain_pt(1), 3)                        +
        4                                           );

    output_pt(0) = retval;
    return;
}

// evaluate f17 function
template<typename T>
void f17(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    if (domain_pt.size() != 3 || output_pt.size() != 1)
    {
        fprintf(stderr, "Error: f17 function is only defined for a 3d domain and scalar output.\n");
        exit(0);
    }

    T E         = domain_pt(0);
    T G         = domain_pt(1);
    T M         = domain_pt(2);
    T gamma     = sqrt(M * M * (M * M + G * G));
    T kprop     = (2.0 * sqrt(2.0) * M * G * gamma ) / (M_PI * sqrt(M * M + gamma));
    T retval    = kprop / ((E * E - M * M) * (E * E - M * M) + M * M * G * G);

    output_pt(0) = retval;
    return;
}

// evaluate f18 function
template<typename T>
void f18(const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs&  args, int k)
{
    if (domain_pt.size() != 4 || output_pt.size() != 1)
    {
        fprintf(stderr, "Error: f18 function is only defined for a 4d domain and scalar output.\n");
        exit(0);
    }

    T x1        = domain_pt(0);
    T x2        = domain_pt(1);
    T x3        = domain_pt(2);
    T x4        = domain_pt(3);
    T retval    = (atanh(x1) + atanh(x2) + atanh(x3) + atanh(x4)) / ((pow(x1, 2) - 1) * pow(x2, -1));

    output_pt(0) = retval;
    return;
}

template<typename T>
void evaluate_function(string fun, const VectorX<T>& domain_pt, VectorX<T>& output_pt, DomainArgs& args, int k)
{
    if (fun == "sine")          return sine(        domain_pt, output_pt, args, k);
    else if (fun == "cosine")   return cosine(      domain_pt, output_pt, args, k);
    else if (fun == "ncosp1")   return ncosp1(      domain_pt, output_pt, args, k);
    else if (fun == "sinc")     return sinc(        domain_pt, output_pt, args, k);
    else if (fun == "psinc1")   return polysinc1(   domain_pt, output_pt, args, k);
    else if (fun == "psinc2")   return polysinc2(   domain_pt, output_pt, args, k);
    else if (fun == "psinc3")   return polysinc3(   domain_pt, output_pt, args, k);
    else if (fun == "psinc4")   return polysinc4(   domain_pt, output_pt, args, k);
    else if (fun == "ml")       return ml(          domain_pt, output_pt, args, k);
    else if (fun == "f16")      return f16(         domain_pt, output_pt, args, k);
    else if (fun == "f17")      return f17(         domain_pt, output_pt, args, k);
    else if (fun == "f18")      return f18(         domain_pt, output_pt, args, k);
    else if (fun == "sinc_sum") return sinc_sum(    domain_pt, output_pt, args);
    else if (fun == "sinc_sum_2") return sinc_sum_2(domain_pt, output_pt, args);
    else if (fun == "ackley")    return ackley(      domain_pt, output_pt, args);
    else if (fun == "schwefel")  return schwefel(   domain_pt, output_pt, args);
    else if (fun == "rastrigin") return rastrigin(  domain_pt, output_pt, args);
    else if (fun == "gaussian_pair1") return gaussian_pair1(  domain_pt, output_pt);
    else if (fun == "gaussian_pair2") return gaussian_pair2(  domain_pt, output_pt);
    else if (fun == "gaussian_mixture") return gaussian_mixture(domain_pt, output_pt);
    else if (fun == "rotating_gaussian") return rotating_gaussian(domain_pt, output_pt);
    else
    {
        cerr << "Invalid function name in evaluate_function. Aborting." << endl;
        exit(0);
    }

    return;
}




#endif // _MFA_EX_FNS
