#include <iostream>
#include <math.h>

# define M_PI           3.14159265358979323846  /* pi */

#include <cassert>
#ifdef NDEBUG
#define assert(condition) ((void)0)
#else
#define assert(condition) /*implementation defined*/
#endif

typedef double RealType;

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

	RealType airfoilY(RealType x); // this is the positive y(x) - upper airfoil - assuming the lower airfoil is (-1)*y(x)
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
	int imax = 121;
	int jmax = 41;
	int iTEL = 21;
	int iLE = 61;
	int iTEU = 101;
	RealType XSF = 1.1;
	RealType YSF = 1.1;
	RealType deltaY = 0.01;
	RealType t = 0.12;
	RealType xint = 1; // 1 for classical 1.008930411365 for modified

	RealType airfoilY(RealType x) {
		return 5 * t * (0.2969 * sqrt(xint * x) - 0.1260 * (xint * x) - 0.3516 * pow(xint * x, 2) + 0.2843 * pow(xint * x, 3) - 0.1015 * pow(xint * x, 4));
	}
};

struct MeshTransform {
	MeshTransformParams *params;
	RealType **S; // 2-D array (RHS to be)
	RealType deltaX;

	// constructor
	MeshTransform(MeshTransformParams* p) {
		S = new RealType *[params->imax];
		for (int i = 0; i < params->imax; i++) {
			S[i] = new RealType[params->jmax];
		}
		deltaX = 1.0 / (params->iLE - params->iTEL); 
	}

	// destructor 
	~MeshTransform() {
		for (int i = 0; i < params->imax; i += 1) {
			delete[] S[i];
		}
		delete[] S;
	}
	
	void initBoundaryCond() {
		for (int i = params->iTEL - 1; i < params->iLE; i++) { // note that C++ indices starts at 0
			x[i][0] = 1.0 - cos(M_PI * (iLE - i) * deltaX / 2.0);
			y[i][0] = airfoilY(x[i][0]);
		}
	}
}