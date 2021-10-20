#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include "Vector3D.h"

using namespace std;
using namespace VECTOR3D;

#define DO_MOL for (int n=0; n<nMol; ++n)

const double L = 20.;
const double rRange = 5.;

double RandomGenerator(double b1, double b2)
{
	const long int a = 16807;
	const long int c = 0;
	const long m = 2147483647;
	static long int x0 = 1;
	long int x1;
	double r;
    
	x1 = (a * x0 + c) % m;
	r = ((double) x1) / ((double) m);
	x0 = x1;
	return b1 + r * (b2 - b1);
}

void PeriodicBC (Vector & r)
{
	double xNew = r.x();
	double yNew = r.y();
	double zNew = r.z();
	if (r.x() >= 0.5 * L) xNew = r.x() - L;
	else if (r.x() < -0.5 * L) xNew = r.x() + L;
	if (r.y() >= 0.5 * L) yNew = r.y() - L;
	else if (r.y() < -0.5 * L) yNew = r.y() + L;
	if (r.z() >= 0.5 * L) zNew = r.z() - L;
	else if (r.z() < -0.5 * L) zNew = r.z() + L;
	r.SetValues(xNew, yNew, zNew);
}

void Initialize (int nMol, double dt, double tempStart, Vector* rO, Vector* rN, Vector* v, Vector* a)
{
	double sigmaMB = sqrt(tempStart);
	double pMax = 1. / (sqrt(2.*M_PI) * sigmaMB);
	double vRand, vTest, p;
	double inf = -5. * sigmaMB;
	double sup = 5. * sigmaMB;
	int molCount = 0;
	int coordCount = 0;
	double coordVal[3];
	
	while (molCount < nMol)
	{
		vRand = RandomGenerator(inf, sup);
		vTest = RandomGenerator(0., 1.);
		p = pMax * exp(-0.5 * pow(vRand / sigmaMB, 2));
		if (vTest < p / pMax)
		{
			coordVal[coordCount] = vRand;
			++ coordCount;
		}
		if (coordCount == 2)
		{
			coordCount = 0;
			v[molCount].SetValues(coordVal[0], coordVal[1], coordVal[2]);
			++ molCount;
		}
	}
	
	double x, y, z;
	DO_MOL
	{
		x = RandomGenerator(-0.5 * L, 0.5 * L);
		y = RandomGenerator(-0.5 * L, 0.5 * L);
		z = RandomGenerator(-0.5 * L, 0.5 * L);
		rO[n].SetValues(x, y, z);
		rN[n] = rO[n] + v[n] * dt;
		PeriodicBC (rN[n]);
		a[n].SetValues(0., 0., 0.);
	}
}

void DoStep (int nMol, double dt, Vector* rO, Vector* rN, Vector* v, Vector* a, double & vvSum)
{
	Vector rCopy;
	vvSum = 0.;
	DO_MOL
	{
		rCopy = rN[n];
		rN[n] = 2. * rCopy - rO[n] + (a[n] * dt * dt);
		PeriodicBC (rN[n]);
		v[n] = (rN[n] - rO[n]) / (2. * dt);
		rO[n] = rN[n];
		vvSum += v[n] * v[n];
	}
}

void ComputeAcc (int nMol, Vector* rN, Vector* a, double & uSum)
{
	Vector dr;
	double r, fVal, rri, rri3;
	DO_MOL a[n].SetValues (0., 0., 0.);
	uSum = 0.;
	for (int i=0; i< nMol-1; ++i)
	{
		for (int j=i+1; j<nMol; ++j)
		{
			dr = rN[i] - rN[j];
			PeriodicBC (dr);
			r = dr.Mag();
			if (r <= rRange)
			{
				rri = pow(1./r,2);
				rri3 = rri * rri * rri;
				fVal = 48. * rri3 * (rri3 - 0.5) * rri;
				a[i] += fVal * dr;
				a[j] += -fVal * dr;
				uSum += 4. * rri3 * (rri3 - 1.);
			}
		}
	}
}

void PrintProp (ostream & out, double timeNow, int nMol, double vvSum, double uSum)
{
	double energy = (uSum + 0.5 * vvSum) / nMol;
	double temperature = vvSum / (3. * nMol);
	out << timeNow << "\t" << temperature << "\t" << energy << endl;
}

void ComputeProp (int nMol, double vvSum, double uSum, double & eSum, double & eeSum, double & tSum, double & ttSum)
{
	double energy = (uSum + 0.5 * vvSum) / nMol;
	double temperature = 1. / (3. * nMol) * vvSum;
	eSum += energy;
	eeSum += energy * energy;
	tSum += temperature;
	ttSum += temperature * temperature;
}

void ComputePCF (int nMol, int* histPCF, Vector* rN, double rMaxPCF, double deltaR)
{
	Vector dr;
	double r;
	int index;
	for (int i=0; i<nMol-1; ++i)
	{
		for (int j=i+1; j<nMol; ++j)
		{
			dr = rN[i] - rN[j];
			PeriodicBC (dr);
			r = dr.Mag();
			if (r < rMaxPCF)
			{
				index = (int) (r / deltaR);
				++ histPCF[index];
			}
		}
	}
}

int main (int argc, char const* argv[])
{
	time_t tStart, tEnd;
	time (& tStart);
	double vol = L * L * L;
	double density, tempStart, dt;
	cout << "Temperatura: ";
	cin >> tempStart;
	cout << "Densita': ";
	cin >> density;
	cout << "dt: ";
	cin >> dt;
	int nMol = density * vol;
	const int stepLimit = 10000;
	const int startAccum = 100;
	
	const int histPCFSize = 50;
	double rMaxPCF = 0.5 * L;
	double deltaR = rMaxPCF / (double) histPCFSize;
	
	Vector* rN = new Vector[nMol];
	Vector* rO = new Vector[nMol];
	Vector* v = new Vector[nMol];
	Vector* a = new Vector[nMol];
	double vvSum = 0., uSum = 0.;
	double eSum = 0., eeSum = 0., tSum = 0., ttSum = 0.;
	double timeNow;
	int* histPCF = new int[histPCFSize];
	
	
	ofstream out1 ("Result.dat");
	out1.precision (10);
	Initialize (nMol, dt, tempStart, rO, rN, v, a);
	int stepCount = 0;
	int accumCount = 0;
	int PCFcount = 0;
	while (stepCount < stepLimit)
	{
		++ stepCount;
		if (stepCount % 1000 == 0)
			cout << "step: " << stepCount << " / " << stepLimit << endl;
		timeNow =  dt * stepCount;
		ComputeAcc (nMol, rN, a, uSum);
		DoStep (nMol, dt, rO, rN, v, a, vvSum);
		if (stepCount > startAccum)
		{
			++ accumCount;
			PrintProp (out1, timeNow, nMol, vvSum, uSum);
			ComputeProp (nMol, vvSum, uSum, eSum, eeSum, tSum, ttSum);
			if (stepCount % 100 == 0)
			{
				++ PCFcount;
				ComputePCF(nMol, histPCF, rN, rMaxPCF, deltaR);
			}
		}
	}
	out1.close();
	
	ofstream out2 ("ResultPCF.dat");
	out2.precision (10);
	double histPCFDouble;
	for (int i=1; i<histPCFSize; ++i)
	{
		histPCFDouble = (double) histPCF[i] * vol / ((double) pow(nMol * i * deltaR,2) * 4.* M_PI * PCFcount);
		out2 << i * deltaR << "\t" << histPCFDouble << endl;
	}
	out2.close();
	
	double eMed, eErr, tMed, tErr;
	eMed = eSum / accumCount;
	eErr = sqrt(((eeSum / accumCount) - pow(eSum / accumCount,2)) / accumCount);
	tMed = tSum / accumCount;
	tErr = sqrt(((ttSum / accumCount) - pow(tSum / accumCount,2)) / accumCount);
	cout << "In Result.dat: ";
	cout << "| 1: time | 2: temperature | 3: EMed |" << endl;
	cout << "Energy: " << eMed << " +- " << eErr << endl;
	cout << "Temperature: " << tMed << " +- " << tErr << endl;
	time (& tEnd);
	cout << "Done: execution time " << (tEnd - tStart) << " s\a" << endl;
	delete [] rN;
	delete [] rO;
	delete [] v;
	delete [] a;
	delete [] histPCF;	
	
	return 0;
}
