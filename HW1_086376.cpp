#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi  under Windows */
#endif

// REAL_TYPE is defined either as float or double (two different programs)
// in order to see impact of percision
typedef float REAL_TYPE;

#define ERROR 1
#define NOERROR 0

/*   a, b, c, are the vectors of the diagonal and the two off-diagonals.
     The vector d is the RHS vector, the vector u is the solution
     vector, "is" is the starting point, and ie is
     the last point.
*/
int tridiag(REAL_TYPE* a, REAL_TYPE* b, REAL_TYPE* c, REAL_TYPE* d, REAL_TYPE* u, int is, int ie)
{
    int i;
    REAL_TYPE beta;

    // ofstream f("tridiag_dump.csv");
    // f << "a[i],b[i],c[i],d[i],beta,new b[i],new d[i]" << endl;
    for (i = is + 1; i <= ie; i++)
    {
        if (b[i - 1] == 0.) return(ERROR);
        beta = a[i] / b[i - 1];
        // f << a[i] << ", " << b[i] << ", " << c[i] << ", " << d[i] << ", " << beta;
        b[i] = b[i] - c[i - 1] * beta;
        d[i] = d[i] - d[i - 1] * beta;
        // f << ", " << b[i] << ", " << d[i] << endl;
    }

    u[ie] = d[ie] / b[ie];
    // f << "u[ie] = " << u[ie] << endl;
    for (i = ie - 1; i >= is; i--)
    {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
        // f << "u[" << i << "] = " << u[i] << endl;
    }
    return(NOERROR);
}

// abstract representation for a second order differential equations as
//    a(x)*y'' + b(x)*y' + c(x)*y = d(x)
struct SecOrderDiffEq {
    virtual REAL_TYPE a(REAL_TYPE x) = 0;
    virtual REAL_TYPE b(REAL_TYPE x) = 0;
    virtual REAL_TYPE c(REAL_TYPE x) = 0;
    virtual REAL_TYPE d(REAL_TYPE x) = 0;
};

// ex0 differntial equation
struct Ex0SecOrderDiffEq : public SecOrderDiffEq {
    REAL_TYPE a(REAL_TYPE x) { return 1.0; }
    REAL_TYPE b(REAL_TYPE x) { return x; }
    REAL_TYPE c(REAL_TYPE x) { return (x * x); }
    REAL_TYPE d(REAL_TYPE x) {
        x *= (2 * M_PI); // M_PI is taken from math.h
        return (sin(x) + cos(x));
    }
};

struct EquiSpaceCellsSpec {
    EquiSpaceCellsSpec(float s_v, float e_v, int n_v) : s(s_v), e(e_v), n(n_v) {}

    float s;   // interval start
    float e;   // interval end
    int n;      // number of equispaced cells
};

// LinearSystem is a base class for Neuman or Dirichlet driven linear systems representaions
// constructed with SecOrderDiffEq object and GidSpec
// it solves the system on construction and provides access to y(i)
class LinearSystem {
public:
    LinearSystem(EquiSpaceCellsSpec& g, SecOrderDiffEq* sode) : gs(g), de(sode) {
        A = (REAL_TYPE*)malloc((gs.n + 1) * sizeof(REAL_TYPE));
        B = (REAL_TYPE*)malloc((gs.n + 1) * sizeof(REAL_TYPE));
        C = (REAL_TYPE*)malloc((gs.n + 1) * sizeof(REAL_TYPE));
        D = (REAL_TYPE*)malloc((gs.n + 1) * sizeof(REAL_TYPE));
        Y = (REAL_TYPE*)malloc((gs.n + 1) * sizeof(REAL_TYPE));

        // calculate the basic A,B,C,D vectors before tuning them for specific method
        h = (gs.e - gs.s) / gs.n;
        for (int i = 0; i <= gs.n; i++) {
            REAL_TYPE x = gs.s + i * h;
            A[i] = (de->a(x) / (h * h)) - (de->b(x) / (2.0 * h));
            B[i] = -(2.0 * de->a(x) / (h * h)) + de->c(x);
            C[i] = de->a(x) / (h * h) + (de->b(x) / (2.0 * h));
            D[i] = de->d(x);
        }
    }
    ~LinearSystem() {
        free(A);
        free(B);
        free(C);
        free(D);
        free(Y);
    }

    virtual void solve() {};   // here tweak the A,B,C,D vectors and call tridiag

    REAL_TYPE y(int i) { return Y[i]; }

protected:
    EquiSpaceCellsSpec gs;
    SecOrderDiffEq* de;
    REAL_TYPE h;
    REAL_TYPE* A;
    REAL_TYPE* B;
    REAL_TYPE* C;
    REAL_TYPE* D;
    REAL_TYPE* Y;
};

class DirichletLinearSystem : public LinearSystem
{
public:
    DirichletLinearSystem(EquiSpaceCellsSpec& g, SecOrderDiffEq* sode, REAL_TYPE y_sv, REAL_TYPE y_ev)
        : LinearSystem(g, sode), y_s(y_sv), y_e(y_ev) {}

    void solve() {
        int n = this->gs.n;

        D[1]     -= (A[1] * y_s);
        D[n - 1] -= (C[n - 1] * y_e);

        if (tridiag(A, B, C, D, Y, 1, n - 1) == ERROR) {
            cerr << "tridiag failed" << endl;
            exit(1);
        }
    }

private:
    REAL_TYPE y_s;
    REAL_TYPE y_e;
};

class NeumannLinearSystem : public LinearSystem
{
public:
    NeumannLinearSystem(EquiSpaceCellsSpec& g, SecOrderDiffEq* sode, REAL_TYPE y1_sv, REAL_TYPE y1_ev)
        : LinearSystem(g, sode), y1_s(y1_sv), y1_e(y1_ev) {}

    void solve() {
        int        n = gs.n;
        D[0] += (2.0 * h * A[0] * y1_s);
        D[n] -= (2.0 * h * C[n] * y1_e);
        C[0] += A[0];
        A[n] += C[n];

        if (tridiag(A, B, C, D, Y, 0, n) == ERROR) {
            cerr << "tridiag failed" << endl;
            exit(1);
        }
    }

private:
    REAL_TYPE y1_s;
    REAL_TYPE y1_e;
};

// process_input
// gets a spec file name with lines of the following format
// <method> <is> <ie> <n> <s_cond> <e_cond> <output-file> <output-step>
// for each line it solves the system for Neuman/Dirichlet with double/flow
// writing to <output-file> a line per Xi where i goes in <output-step>
// <output-step> is required for very large n grids where for display we only
// want to a relatively small sample of points (e.g. 500)
void process_input(const char* spec_file_name)
{
    ifstream sf(spec_file_name);
    REAL_TYPE s, e;                     // function support on [s,e]
    REAL_TYPE s_cond, e_cond;           // bounday conditions on s,e
    int n, step;
    string out_file_name;
    string method;

    Ex0SecOrderDiffEq de;

    while (!sf.eof()) {
        sf >> method>> s >> e >> n >> s_cond >> e_cond >> out_file_name >> step;
        if (sf.fail()) {
            break;
        }
        EquiSpaceCellsSpec gs(s, e, n);
        REAL_TYPE h = (gs.e - gs.s) / gs.n;

        ofstream of(out_file_name);
        if (method == "Neumann") {
            NeumannLinearSystem neumann(gs, &de, s_cond, e_cond);
            neumann.solve();

            of << "Neumann on [" << gs.s << ", " << gs.e << "] with Y'(" << gs.s << ")=" << s_cond
                << " and Y'(" << gs.e << ")=" << e_cond << endl;

            of << "X, Neumann" << endl;
            for (int i = 0; i <= gs.n; i += step) {
                of << (gs.s + i * h) << ", " << neumann.y(i) << endl;
            }
        }
        else if (method == "Dirichlet") {
            DirichletLinearSystem dirichlet(gs, &de, s_cond, e_cond);
            dirichlet.solve();

            of << "Dirichlet on [" << gs.s << ", " << gs.e << "] with Y(" << gs.s << ")=" << s_cond
                << " and Y(" << gs.e << ")=" << e_cond << endl;

            of << "X, Dirichlet" << endl;
            of << gs.s << ", " << s_cond << endl;
            for (int i = step; i < gs.n; i += step) {
                of << (gs.s + i * h) << ", " << dirichlet.y(i) << endl;
            }
            of << gs.e << ", " << e_cond << endl;
        }
        else {
            cerr << "unknown solving method! chiaoooooooooooo" << endl;
            exit(1);
        }
    }
}

int main(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i) {
        process_input(argv[i]);
    }
    return 0;
}
