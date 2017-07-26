/* Numerical solvers for quadratic, cubic, and quartic equations.
 * All solvers return the number of real roots.
 */

#include <math.h>

// y = sqrt(x)
void sqrt_cplx(double y[2], double x[2]) {
    double m = pow(x[0]*x[0]+x[1]*x[1], 0.25);
    double th = 0.5*atan2(x[1], x[0]);
    y[0] = m*cos(th);
    y[1] = m*sin(th);
}

// a x^2 + b x + c (following NR in C)
int solve_quadratic(double x[2], double a, double b, double c) {
    double d = b*b-4*a*c;
    if(d < 0.0) { // cplx
        double q0 = -0.5*b;
        double q1 = 0.5*sqrt(-d);
        double q2 = q0*q0+q1*q1;
        if(a*a > q2) { // test which to divide by
            x[0] = q0/a;
            x[1] = q1/a;
        } else {
            x[0] = q0*c/q2;
            x[1] = q1*c/q2;
        }
        return 0;
    } else {
        double q;
        if(b < 0.0) {
           q = -0.5*(b - sqrt(d));
        } else {
           q = -0.5*(b + sqrt(d));
        }
        x[0] = q/a;
        x[1] = c/q;
        return 2;
    }
}

int quadratic(sil_State *S) {
    double x[2];
    double b = sil_todouble(S, 1);
    double c = sil_todouble(S, 2);
    int n = solve_quadratic(x, 1.0, b, c);
    sil_settop(S, 0);
    sil_pushinteger(S, n);
    sil_pushdouble(S, x[0]);
    sil_pushdouble(S, x[1]);
    sil_settuple(S, 3);
    return 0;
}

// x^3 + a x^2 + b x + c = 0
int solve_cubic(double x[3], double a, double b, double c) {
    if(c == 0.0) { // handle degenerate quadratic
        x[0] = 0.0;
        return 1+solve_quadratic(x+1, 1.0, a, b);
    }
    double Q = (a*a - 3*b)/9.0;
    double R = (2*a*a*a - 9*a*b + 27*c)/54.0;
    double d = R*R - Q*Q*Q;

    if(d < 0.0) { // 3 real roots.
        double th = acos(R/sqrt(Q*Q*Q));
        Q = -2*sqrt(Q);
        x[0] = Q*cos(th/3.0) - a/3.0;
        x[1] = Q*cos((th + 2*M_PI)/3.0) - a/3.0;
        x[2] = Q*cos((th - 2*M_PI)/3.0) - a/3.0;
        return 3;
    }
    double A = cbrt(fabs(R) + sqrt(d));
    if(R >= 0.0) A = -A; // * -sgn(R)
    double B = 0.0;
    if(A != 0.0) {
        B = Q/A;
    }
    x[0] = A + B - a/3.0;
    x[1] = -0.5*(A+B) - a/3.0;  // real part of 2nd root
    x[2] = 0.5*sqrt(3.0)*(A-B); // cplx part
    return 1;
}

int cubic(sil_State *S) {
    double x[3];
    double a = sil_todouble(S, 1);
    double b = sil_todouble(S, 2);
    double c = sil_todouble(S, 3);
    int n = solve_cubic(x, a, b, c);
    sil_settop(S, 0);
    sil_pushinteger(S, n);
    sil_pushdouble(S, x[0]);
    sil_pushdouble(S, x[1]);
    sil_pushdouble(S, x[2]);
    sil_settuple(S, 4);
    return 0;
}

// Solve:
// y^4 + p y^2 + q y + r = 0
// By factoring into (y^2 + sy + t) (y^2 + uy + v)
int solve_depressed(double x[4], double p, double q, double r) {
    double z[3];
    int n;

    if(r == 0.0) { // y = 0 is a trivial root
        x[0] = 0.0;
        return solve_cubic(x+1, 0.0, p, q) + 1;
    }
    if(fabs(q) < 1e-15) { // biquadratic eqn.
biquad:
        n = solve_quadratic(z, 1.0, p, r);
        if(n == 0) {
            sqrt_cplx(x, z);
            x[2] = -x[0];
            x[3] =  x[1];
            return 0;
        }
        if(z[1] > z[0]) { // space by 2, put largest first
            x[0] = z[1];
            x[2] = z[0];
        } else {
            x[0] = z[0];
            x[2] = z[1];
        }
        n = 0;
        if(x[0] >= 0.0) {
            x[0] = sqrt(x[0]);
            x[1] = -x[0];
            n += 2;
        } else {
            x[1] = 0.0;
            sqrt_cplx(x, x);
        }
        if(x[2] >= 0.0) {
            x[2] = sqrt(x[2]);
            x[3] = -x[2];
            n += 2;
        } else {
            x[3] = 0.0;
            sqrt_cplx(x+2, x+2);
        }
        return n;
    }
    n = solve_cubic(z, 2*p, p*p-4*r, -q*q);
    //printf("resolvent, %d: %e %e %e\n", n, z[0], z[1], z[2]);
    if(n == 3) { // Use the largest root.
        z[0] = z[1] > z[0] ? z[1] : z[0];
        z[0] = z[2] > z[0] ? z[2] : z[0];
    }
    if(z[0] < 0.0) {
        if(fabs(z[0]) < 1e-12) { // should have solved bi-quadratic.
            goto biquad;
        }
        printf("Unable to find a root: %e %e %e, %e\n", p, q, r, z[0]);
        return 0;
    }
    double u = sqrt(z[0]);
    double s = -u;
    double t = 0.5*(p + u*u + q/u);
    double v = 0.5*(p + u*u - q/u);
    //printf("(x^2 + %e x + %e) (x^2 + %e x + %e)\n", s, t, u, v);
    n = solve_quadratic(x, 1.0, s, t);
    if(n == 0) {
        x[2] = x[0];
        x[3] = x[1];
        return solve_quadratic(x, 1.0, u, v);
    }
    return n + solve_quadratic(x+2, 1.0, u, v);
}

// a x^4 + b x^3 + c x^2 + d x + e = 0
// (x = y - b/4a)
int solve_quartic(double x[4], double a, double b, double c,
                               double d, double e) {
    double p = (8*a*c - 3*b*b)/(8*a*a);
    double q = (b*b*b - 4*a*b*c + 8*a*a*d)/(8*a*a*a);
    double r = (64*a*a*(4*a*e - b*d) + b*b*(16*c*a-3*b*b))/(256*a*a*a*a);

    //printf("p=%e q=%e r=%e\n", p, q, r);
    int i, n = solve_depressed(x, p, q, r);

    // Shift solutions
    double x0 = -0.25*b/a;
    for(i=0; i<n; i++) {
        x[i] += x0;
    }
    for(; i<4; i+=2) {
        x[i] += x0;
    }
    return n;
}

int quartic(sil_State *S) {
    double x[4];
    double b = sil_todouble(S, 1);
    double c = sil_todouble(S, 2);
    double d = sil_todouble(S, 3);
    double e = sil_todouble(S, 4);
    int n = solve_quartic(x[4], 1.0, b, c, d, e);

    sil_settop(S, 0);
    sil_pushinteger(S, n);
    sil_pushdouble(S, x[0]);
    sil_pushdouble(S, x[1]);
    sil_pushdouble(S, x[2]);
    sil_pushdouble(S, x[3]);
    sil_settuple(S, 5);
    return 0;
}

