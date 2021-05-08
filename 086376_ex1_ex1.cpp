#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm> 

# define M_PI           3.14159265358979323846  /* pi */

#include <cassert>

typedef double RealType;

int btri2s(RealType* a, RealType* b, RealType* c, RealType* f, int kd, int ks, int ke)
{
	/* Local variables */
	int k, kr;
	int md = 2;
	RealType b11, b12, b21, b22, beta11, beta12, beta21, beta22;
	RealType det, rdet, f1, f2;

	/*   (A,B,C)F = F, F and B are overloaded, solution in F */

	/*   Part 1. Forward block sweep */

	for (k = ks + 1; k <= ke; k++)
	{

		/* Invert B */

		kr = k - 1;
		det = b[kr + kd * (0 + md * 0)] * b[kr + kd * (1 + md * 1)] \
			- b[kr + kd * (0 + md * 1)] * b[kr + kd * (1 + md * 0)];
		if (det == 0.0) return (-1);
		rdet = 1.0 / det;
		b11 = b[kr + kd * (1 + md * 1)] * rdet;
		b12 = -b[kr + kd * (0 + md * 1)] * rdet;
		b21 = -b[kr + kd * (1 + md * 0)] * rdet;
		b22 = b[kr + kd * (0 + md * 0)] * rdet;

		/* beta = a(k) * binv(k-1) */

		beta11 = a[k + kd * (0 + md * 0)] * b11 + a[k + kd * (0 + md * 1)] * b21;
		beta12 = a[k + kd * (0 + md * 0)] * b12 + a[k + kd * (0 + md * 1)] * b22;
		beta21 = a[k + kd * (1 + md * 0)] * b11 + a[k + kd * (1 + md * 1)] * b21;
		beta22 = a[k + kd * (1 + md * 0)] * b12 + a[k + kd * (1 + md * 1)] * b22;

		/* b(k) = b(k) - beta * c(k-1) */

		b[k + kd * (0 + md * 0)] = b[k + kd * (0 + md * 0)] \
			- beta11 * c[k - 1 + kd * (0 + md * 0)] \
			- beta12 * c[k - 1 + kd * (1 + md * 0)];
		b[k + kd * (0 + md * 1)] = b[k + kd * (0 + md * 1)] \
			- beta11 * c[k - 1 + kd * (0 + md * 1)] \
			- beta12 * c[k - 1 + kd * (1 + md * 1)];
		b[k + kd * (1 + md * 0)] = b[k + kd * (1 + md * 0)] \
			- beta21 * c[k - 1 + kd * (0 + md * 0)] \
			- beta22 * c[k - 1 + kd * (1 + md * 0)];
		b[k + kd * (1 + md * 1)] = b[k + kd * (1 + md * 1)] \
			- beta21 * c[k - 1 + kd * (0 + md * 1)] \
			- beta22 * c[k - 1 + kd * (1 + md * 1)];

		/* f(k) = f(k) - beta * f(k-1) */

		f[k + kd * 0] = f[k + kd * 0] \
			- beta11 * f[k - 1 + kd * 0] \
			- beta12 * f[k - 1 + kd * 1];
		f[k + kd * 1] = f[k + kd * 1] \
			- beta21 * f[k - 1 + kd * 0] \
			- beta22 * f[k - 1 + kd * 1];

	}

	/* Backward substitution */

	k = ke;
	det = b[k + kd * (0 + md * 0)] * b[k + kd * (1 + md * 1)] \
		- b[k + kd * (0 + md * 1)] * b[k + kd * (1 + md * 0)];
	if (det == 0.0) return (-1);
	rdet = 1.0 / det;
	b11 = b[k + kd * (1 + md * 1)] * rdet;
	b12 = -b[k + kd * (0 + md * 1)] * rdet;
	b21 = -b[k + kd * (1 + md * 0)] * rdet;
	b22 = b[k + kd * (0 + md * 0)] * rdet;

	f1 = b11 * f[k + kd * 0] + b12 * f[k + kd * 1];
	f2 = b21 * f[k + kd * 0] + b22 * f[k + kd * 1];

	f[k + kd * 0] = f1;
	f[k + kd * 1] = f2;

	for (k = ke - 1; k >= ks; --k)
	{

		det = b[k + kd * (0 + md * 0)] * b[k + kd * (1 + md * 1)] \
			- b[k + kd * (0 + md * 1)] * b[k + kd * (1 + md * 0)];
		if (det == 0.0) return (-1);
		rdet = 1.0 / det;
		b11 = b[k + kd * (1 + md * 1)] * rdet;
		b12 = -b[k + kd * (0 + md * 1)] * rdet;
		b21 = -b[k + kd * (1 + md * 0)] * rdet;
		b22 = b[k + kd * (0 + md * 0)] * rdet;

		f[k + kd * 0] = f[k + kd * 0] \
			- c[k + kd * (0 + md * 0)] * f[k + 1 + kd * 0] \
			- c[k + kd * (0 + md * 1)] * f[k + 1 + kd * 1];
		f[k + kd * 1] = f[k + kd * 1] \
			- c[k + kd * (1 + md * 0)] * f[k + 1 + kd * 0] \
			- c[k + kd * (1 + md * 1)] * f[k + 1 + kd * 1];

		f1 = b11 * f[k + kd * 0] + b12 * f[k + kd * 1];
		f2 = b21 * f[k + kd * 0] + b22 * f[k + kd * 1];

		f[k + kd * 0] = f1;
		f[k + kd * 1] = f2;
	}

	return 0;

} /* btri4s_ */

void debugPrintMat(const char *outFileName, RealType** m, int rows, int cols) {
	std::ofstream out(outFileName);
	for (int j = rows - 1; j >= 0; j--) {
		for (int i = 0; i < cols; i++) {
			out << m[i][j] << ", ";
		}
		out << std::endl;
	}
}

// the general definition of the problems parameters:
struct MeshTransformParams {
	int imax;
	int jmax;
	int iTEL;
	int iLE;
	int iTEU;
	RealType XSF;
	RealType YSF;
	RealType deltaY;
	RealType t;
	RealType xint;
	RealType deltaEta;
	RealType deltaPsi;

	virtual RealType airfoilY(RealType x) = 0; // this is the positive y(x) - upper airfoil - assuming the lower airfoil is (-1)*y(x)
	void checkParams() {
		assert(iTEL > 0);
		assert(iTEU > 0);
		assert(iLE > iTEL);
		assert(jmax > 0);
		assert(imax > 0);
	}
};

// specific parameters for this exrecise :)
struct Ex1MeshTransformParam : public MeshTransformParams {
	Ex1MeshTransformParam() {
		imax = 121;
		jmax = 41;
		iTEL = 21;
		iLE = 61;
		iTEU = 101;
		XSF = 1.1;
		YSF = 1.1;
		deltaY = 0.01;
		t = 0.12;
		xint = 1; // 1 for classical 1.008930411365 for modified
		deltaEta = 1;
		deltaPsi = 1;
	}

	virtual RealType airfoilY(RealType x) {
		return 5 * t * (0.2969 * sqrt(xint * x) - 0.1260 * (xint * x) - 0.3516 * pow(xint * x, 2) + 0.2843 * pow(xint * x, 3) - 0.1015 * pow(xint * x, 4));
	}
};

struct MeshTransform {
	MeshTransformParams* params;
	RealType** arr; // 2-D array (RHS to be)
	RealType** x;
	RealType** y;
	RealType** xEta;
	RealType** xPsi;
	RealType** yEta;
	RealType** yPsi;
	RealType deltaX;
	RealType** deltaS;
	RealType A[2][2];
	RealType B[2][2];
	RealType S[2];
	RealType Z[2];
	RealType* a;
	RealType* b;
	RealType* c;
	RealType* f;

	// constructor
	MeshTransform(MeshTransformParams* p) : params(p) {
		arr = new RealType * [params->imax];
		x = new RealType * [params->imax];
		y = new RealType * [params->imax];
		xEta = new RealType * [params->imax];
		xPsi = new RealType * [params->imax];
		yEta = new RealType * [params->imax];
		yPsi = new RealType * [params->imax];
		deltaS = new RealType * [params->imax];
		a = new RealType[4 * params->imax];
		b = new RealType[4 * params->imax];
		c = new RealType[4 * params->imax];
		f = new RealType[2 * params->imax];

		for (int i = 0; i < params->imax; i++) {
			arr[i] = new RealType[params->jmax];
			x[i] = new RealType[params->jmax];
			y[i] = new RealType[params->jmax];
			xEta[i] = new RealType[params->jmax];
			xPsi[i] = new RealType[params->jmax];
			yEta[i] = new RealType[params->jmax];
			yPsi[i] = new RealType[params->jmax];
			deltaS[i] = new RealType[params->jmax];
			std::fill_n(x[i], params->jmax, 0);
		}
		deltaX = 1.0 / (params->iLE - params->iTEL);
	}

	// destructor 
	~MeshTransform() {
		for (int i = 0; i < params->imax; i++) {
			delete[] arr[i];
			delete[] x[i];
			delete[] y[i];
			delete[] xEta[i];
			delete[] xPsi[i];
			delete[] yEta[i];
			delete[] yPsi[i];
			delete[] deltaS[i];
		}
		delete[] arr;
		delete[] x;
		delete[] y;
		delete[] xEta;
		delete[] xPsi;
		delete[] yEta;
		delete[] yPsi;
		delete[] deltaS;
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] f;
	}

	void initBoundaryCond() {
		// all defined in the guidance file...
		// to seperate to little functions...
		// lower + upper (7)
		for (int i = params->iTEL - 1; i < params->iLE; i++) { // note that C++ indices starts at 0
			x[i][0] = 1.0 - cos(M_PI * RealType(params->iLE - i) * deltaX / 2.0);
			y[i][0] = params->airfoilY(x[i][0]) * (-1);
		}
		for (int i = params->iLE - 1; i < params->iTEU; i++) { // note that C++ indices starts at 0
			x[i][0] = 1.0 - cos(M_PI * RealType(params->iLE - i) * deltaX / 2.0);
			y[i][0] = params->airfoilY(x[i][0]);
		}
		// the rest
		// 8
		for (int i = params->iTEU - 1; i < params->imax; i++) {
			x[i][0] = x[i - 1][0] + (x[i - 1][0] - x[i - 2][0]) * params->XSF;
			y[i][0] = 0;
		}
		for (int i = 0; i < params->iTEL; i++) {
			x[i][0] = x[params->imax - i - 1][0];
			y[i][0] = 0;
		}
		// 9
		y[params->imax - 1][1] = params->deltaY;

		for (int j = 2; j < params->jmax; j++) {
			y[params->imax - 1][j] = y[params->imax - 1][j - 1] + (y[params->imax - 1][j - 1] - y[params->imax - 1][j - 2]) * params->YSF;
		}
		for (int j = 1; j < params->jmax; j++) {
			x[params->imax - 1][j] = x[params->imax - 1][0];
		}
		for (int j = 0; j < params->jmax; j++) {
			y[0][j] = -y[params->imax - 1][j];
			x[0][j] = x[0][0];
		}
	}

	void calcPartialDeriv(int j) {
		// 2nd order central difference operator (for xPsi xEta yPsi yEta)
		for (int i = 0; i < params->imax; i++) {
			xPsi[i][j] = 0.5 * (x[i + 1][j] - x[i - 1][j]);
			yPsi[i][j] = 0.5 * (y[i + 1][j] - y[i - 1][j]);
			xEta[i][j] = 0.5 * (x[i][j + 1] - x[i][j - 1]);
			yPsi[i][j] = 0.5 * (y[i][j + 1] - y[i][j - 1]);
			deltaS[i][j] = sqrt(pow(xPsi[i][j], 2) + pow(yPsi[i][j], 2)); // cell area eval	
		}
	}
	void calcLHS_A(int i,int j) {
		A[0][0] = xEta[i][j];
		A[0][1] = yEta[i][j];
		A[1][0] = yEta[i][j];
		A[1][1] = -xEta[i][j];
	}
	void calcLHS_B(int i,int j) {
		B[0][0] = xPsi[i][j];
		B[0][1] = yPsi[i][j];
		B[1][0] = -yPsi[i][j];
		B[1][1] = xPsi[i][j];
	}
	void calcLHS_a(int i, int j) {
		a[j + 0 * params->imax] = -0.5 * A[0][0];
		a[j + 1 * params->imax] = -0.5 * A[0][1];
		a[j + 2 * params->imax] = -0.5 * A[1][0];
		a[j + 3 * params->imax] = -0.5 * A[1][1];
	}
	void calcLHS_b(int i, int j) {
		b[j + 0 * params->imax] = B[0][0];
		b[j + 1 * params->imax] = B[0][1];
		b[j + 2 * params->imax] = B[1][0];
		b[j + 3 * params->imax] = B[1][1];
	}
	void calcLHS_c(int i, int j) {
		c[j + 0 * params->imax] = 0.5 * A[0][0];
		c[j + 1 * params->imax] = 0.5 * A[0][1];
		c[j + 2 * params->imax] = 0.5 * A[1][0];
		c[j + 3 * params->imax] = 0.5 * A[1][1];
	}
	void calcLHS(int j) {
		for (int i = 0; i < params->imax; i++) {
			calcLHS_A(i, j);
			calcLHS_B(i, j);
			//def a,b,c for btri2s
			calcLHS_a(i, j);
			calcLHS_b(i, j);
			calcLHS_c(i, j);
		}
	}
	void calcRHS_S(int i, int j) {
		S[0] = 0;
		S[1] = 0 + (xPsi[i][j] * yEta[i][j] - yPsi[i][j] * xEta[i][j]);
		if (i == 0 || i == params->imax - 1) {
			S[1] = 0.5 * S[1];
		}
	}
	void calcRHS_Z(int i, int j) {
		Z[0] = x[i][j];
		Z[1] = y[i][j];
	}
	void calcRHS_f(int i, int j) {
		// f = S[i][j] + B[i][j] * Z[i][j]
		f[i + 0 * params->imax] = S[0] + B[0][0] * Z[0] + B[0][1] * Z[1];
		f[i + 1 * params->imax] = S[1] + B[1][0] * Z[0] + B[1][1] * Z[1];
	}
	void calcRHS(int j) {
		for (int i = 0; i < params->imax; i++) {
			calcRHS_S(i, j);
			calcLHS_B(i, j);
			calcRHS_Z(i, j);
			calcRHS_f(i, j);
		}
	}
};

int main(int argc, char *argv[]) {
	Ex1MeshTransformParam params;
	MeshTransform mesh(&params);
	mesh.initBoundaryCond();
	debugPrintMat("x.csv", mesh.x, params.jmax, params.imax);
	return 0;
}

