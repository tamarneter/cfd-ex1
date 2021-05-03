#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm> 

# define M_PI           3.14159265358979323846  /* pi */

#include <cassert>

typedef double RealType;

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
	RealType** A;
	RealType** B;
	RealType** S;
	RealType** Z;

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

	void calcPartialDeriv() {
		// 2nd order central difference operator (for xPsi xEta yPsi yEta)
		for (int i = 0; i < params->imax; i++) {
			for (int j = 0; j < params->jmax; j++) {
				xPsi[i][j] = 0.5 * (x[i + 1][j] - x[i - 1][j]);
				yPsi[i][j] = 0.5 * (y[i + 1][j] - y[i - 1][j]);
				xEta[i][j] = 0.5 * (x[i][j + 1] - x[i][j - 1]);
				yPsi[i][j] = 0.5 * (y[i][j + 1] - y[i][j - 1]);
				deltaS[i][j] = sqrt(pow(xPsi[i][j], 2) + pow(yPsi[i][j], 2)); // cell area eval
 			}	
		}
	}

	void calcMarchingScheme() {
		for (int i = 0; i < params->imax; i++) {
			for (int j = 0; j < params->jmax; j++) {
				// def A
				A[1][1] = x[i][j];
				A[1][2] = y[i][j];
				A[2][1] = y[i][j];
				A[2][2] = -x[i][j];
 			}
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

