/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Pbc.h"
#include "../tools/Matrix.h"
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;

#define PI 3.14159265

using namespace std;

namespace PLMD {
namespace colvar {

typedef Eigen::Matrix<dual, 3, 1> Vector;
Vector deltaBend2(Vector v1, Vector v2) {return v2-v1;}
Vector crossBend2(Vector a, Vector b) {
	Vector c(3); 
	c << a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0); 
	return c;
}
dual dot2(Vector a, Vector b){
	dual c = a(0)*b(0)+a(1)*b(1)+a(2)*b(2);
	return c;
}

typedef Eigen::Matrix<dual, 3, 3> EigenMatrix;

//+PLUMEDOC COLVAR BEND
/*
 * BEND colvar is a dsDNA/dsRNA restraint that keeps a certain bend angle
 * between the chosen base pairs in a nucleic acid fragment, and uses an automatic 
 * deferentiation algorythm to calculate the atoms derivatives.
 * As input provide 60 atoms from bases that are restrained (4*4 bases * 3 atoms
 * plus 6*2 auxiliary atoms); a force constant, and desired value of a bend angle.
 * The energy penalty will be added to the potential energy functional is
 * E_bend=0.5*k*(bend0-bend)^2
 */
//+ENDPLUMEDOC

class Bend4bp : public Colvar {
	bool pbc;
	static Vector rotVector6Atoms(int arr[6],Eigen::Matrix <dual, 3, 60> atomcoords);
	static Vector rotVectorBasis(EigenMatrix R1,  EigenMatrix R2);
	static dual rotAngleBasis(EigenMatrix P1,  EigenMatrix P2);
	static dual traceBasis(EigenMatrix P1, EigenMatrix P2);
	static EigenMatrix basis(Vector a, Vector b, Vector c);
	static EigenMatrix basis2Vectors(Vector e1, Vector diff);
	static EigenMatrix basis2Vectors_2(Vector e1, Vector diff);
	static dual bendAngle(Eigen::Matrix <dual, 180, 1> atomcoods1D);
	
public:
	explicit Bend4bp(const ActionOptions&);
	virtual void calculate();
	static void registerKeywords(Keywords& keys);
	Vector getAtomPosition_forBend(int i) {
		PLMD::Vector v = getPosition(i);
		Vector AtomCoordinates(v[0],v[1],v[2]);
		return AtomCoordinates;
	};
};

PLUMED_REGISTER_ACTION(Bend4bp,"BENDDNA4BP")

void Bend4bp::registerKeywords(Keywords& keys) {
	Colvar::registerKeywords(keys);
	keys.add("atoms", "ATOMS", "the keyword with which you specify what atoms to use");
}

Bend4bp::Bend4bp(const ActionOptions&ao):
	PLUMED_COLVAR_INIT(ao),
	pbc(true)
{
	vector<AtomNumber> atoms;
	parseAtomList("ATOMS", atoms);
	if (atoms.size() != 60)
	{
		error("Number of specified atoms should be 60");
		log.printf("number or atoms %d \n ", atoms.size());
	}

	bool nopbc = !pbc;
	parseFlag("NOPBC", nopbc);
	pbc = !nopbc;
	if (pbc) log.printf("  using periodic boundary conditions\n");
	else    log.printf("  without periodic boundary conditions\n");

	addValueWithDerivatives();
	setNotPeriodic();
	requestAtoms(atoms);
	checkRead();
}

// calculator
void Bend4bp::calculate() {

	if (pbc) makeWhole();

	Eigen::Matrix<dual, 3, 60 > atomcoords;
	for (int i = 0; i< 60; i++){
		atomcoords.col(i) = getAtomPosition_forBend(i);
	}

	Eigen::Matrix<dual, 180, 1> atomcoods1D;
	for (int c = 0; c < 60; c++) {
		for (int r = 0; r < 3; r++) {
			atomcoods1D[r+c*3] = atomcoords(r,c);
		}
	}

	dual bend = bendAngle(atomcoods1D);
	dual u;

	Eigen::Matrix<double,180,1> atomsDerivatives1D = gradient(bendAngle,wrt(atomcoods1D), at(atomcoods1D),u);

	Eigen::Matrix<dual,3,60> atomsDerivatives;
	for (int c = 0; c < 60; c++) {
		for (int r = 0; r < 3; r++) {
			atomsDerivatives(r,c) = atomsDerivatives1D[r+c*3];
		}
	}

	for (int atom = 0; atom < 60; atom++) {
		PLMD::Vector derivatives;
		for (int i = 0; i < 3; i++) {
			derivatives(i) -= val(atomsDerivatives(atom,i));
		}
		//printf("Atom %d, derx %g, dery %g, derz %g \n", atom, derivatives(0),derivatives(1),derivatives(2));
		setAtomsDerivatives(atom, derivatives);
	}

	//printf("1st atom derivatives: %g, %g, %g \n", atomDerivatives[0][0], atomDerivatives[0][1],atomDerivatives[0][2]);

	setValue(val(bend));
	setBoxDerivativesNoPbc();
}

dual Bend4bp:: bendAngle(Eigen::Matrix <dual, 180, 1> atomcoods1D)
{
	Eigen::Matrix<dual,3,60> atomcoords;
	for (int c = 0; c < 60; c++) {
		for (int r = 0; r < 3; r++) {
			atomcoords(r,c) = atomcoods1D[r+c*3];
		}
	}

	int atomNmbs[14][6] =
	{
		{ 0, 1, 2, 3, 4, 5},      //Basis 1 - strand 1 begins
		{ 3, 4, 5, 6, 7, 8},      //Basis 2
		{ 6, 7, 8, 9, 10, 11},    //Basis 3
		{12, 13, 14, 15, 16, 17}, //Middle atoms strand #1
		{18, 19, 20, 21, 22, 23}, //Basis 4
		{21, 22, 23, 24, 25, 26}, //Basis 5 
		{24, 25, 26, 27, 28, 29}, //Basis 6
		{30, 31, 32, 33, 34, 35}, //Basis 7 - strand 2 begins 
		{33, 34, 35, 36, 37, 38}, //Basis 8
		{36, 37, 38, 39, 40, 41}, //Basis 9
		{42, 43, 44, 45, 46, 47}, //Middle atoms strand #2
		{48, 49, 50, 51, 52, 53}, //Basis 10
		{51, 52, 53, 54, 55, 56}, //Basis 11
		{54, 55, 56, 57, 58, 59}  //Basis 12
	};

	int Direction[14] = {1, 1, 1, 0, 1, 1, 1, -1, -1, -1, 0, -1, -1, -1};
	int Top_Btm[14]  = {1, 1, 1, 0, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1};

// defining rotation vectors for each of the nucleotides pair

	Vector U; Vector U_avTop; Vector U_avBtm;

	for (int n_basis = 0; n_basis < 14; n_basis++) {
		if (Direction[n_basis] == 1 || Direction[n_basis] == -1) {
			U = rotVector6Atoms(atomNmbs[n_basis],atomcoords);
			U *= 1.0 / U.norm();

			// calculating averaged rotation vectors for top and bottom base pairs
			// also named "handles" in 
			// Magnitude and direction of DNA bending induced by screw-axis orientation: influence of sequence, mismatches and abasic sites
         // Jeremy Curuksu,  Krystyna Zakrzewska,  Martin Zacharias
         // Nucleic Acids Research, Volume 36, Issue 7, 1 April 2008, Pages 2268–2283,
			
			if     (Top_Btm[n_basis] == 1) {
				U_avTop += Direction[n_basis] * U;
			}
			else if (Top_Btm[n_basis] == -1) {
				U_avBtm += Direction[n_basis] * U;
			}

		} // finish if loop
	} // finish loop n_basis

	dual bend = acos(dot2(U_avTop,U_avBtm)/(U_avTop.norm()*U_avBtm.norm()))*180/PI;
	
   //printf("Bend: %g \n ",val(bend));
	return bend;
}

EigenMatrix Bend4bp::basis(Vector a, Vector b, Vector c) {
	return basis2Vectors(deltaBend2(a, b), deltaBend2(a, c));
}

EigenMatrix Bend4bp::basis2Vectors(Vector e1, Vector diff) {
	Vector e2 = crossBend2(e1,diff);
	Vector e3 = crossBend2(e1,e2);
	EigenMatrix M(3, 3);
	for (int i = 0; i < 3; i++)
	{
		M(0, i) = e1(i) / e1.norm();
		M(1, i) = e2(i) / e2.norm();
		M(2, i) = e3(i) / e3.norm();
	}
	return M;
}

// basis2Vectors_2 is a different version of the basis function from 2 vectors
EigenMatrix Bend4bp::basis2Vectors_2(Vector RotVect, Vector N1N9) {
	Vector e1 = crossBend2(RotVect,N1N9);
	Vector e2 = crossBend2(RotVect,e1);
	EigenMatrix M(3, 3);
	for (int i = 0; i < 3; i++)
	{
		M(0, i) = e1(i) / e1.norm();
		M(1, i) = e2(i) / e2.norm();
		M(2, i) = RotVect(i) / RotVect.norm();
	}
	return M;
}

Vector Bend4bp::rotVector6Atoms(int arr[6],Eigen::Matrix <dual, 3, 60> atomcoords) {
	int i0 = arr[0]; int i1 = arr[1]; int i2 = arr[2];
	int i3 = arr[3]; int i4 = arr[4]; int i5 = arr[5];
	EigenMatrix B1 = basis(atomcoords.col(i0), atomcoords.col(i1), atomcoords.col(i2));
	EigenMatrix B2 = basis(atomcoords.col(i3), atomcoords.col(i4), atomcoords.col(i5));
	return rotVectorBasis(B1, B2);
}

Vector Bend4bp::rotVectorBasis( EigenMatrix R1,  EigenMatrix R2) {
	EigenMatrix Q(3, 3);
	EigenMatrix R1T(3, 3);
	R1T = R1.transpose();
	Q = R1T*R2;
	//mult(R1T, R2, Q);
	Vector u;
	u(0) = Q(1, 2) - Q(2, 1);
	u(1) = Q(2, 0) - Q(0, 2);
	u(2) = Q(0, 1) - Q(1, 0);
	//u *= 1/u.norm();
	return u;
}

dual Bend4bp::rotAngleBasis( EigenMatrix P1, EigenMatrix P2) {
	dual trace = traceBasis(P1, P2);
	dual tr = (trace - 1) * 0.5;
	if (tr > 1.0) tr = 1.0;
	else if (tr < -1.0) tr = -1.0;
	dual theta = acos(tr); //*180.0/PI;
	if (theta < 0.0) { fprintf(stderr, "negative theta \n");}
	return theta;
}

dual Bend4bp:: traceBasis( EigenMatrix P1, EigenMatrix P2) {
	EigenMatrix Q(3, 3);
	EigenMatrix P1T(3, 3);
	P1T = P1.transpose();
	Q = P1T*P2;
	dual trace = Q(0, 0) + Q(1, 1) + Q(2, 2);
	return trace;
}

}
}
