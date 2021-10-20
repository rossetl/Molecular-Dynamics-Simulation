#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include "Vector3D.h"

using namespace std;
using namespace VECTOR3D;

#define DO_MOL for (int n=0; n<nMol; ++n)
#define RAND_R RandomGenerator (-rRange, rRange)
#define RAND_V RandomGenerator (-vRange, vRange)

int acceptCount = 0;
int totalCount = 0;
const double rMax = 5.;
double rRange, vRange;
const double L = 10.;

double RandomGenerator(double b1, double b2)
{
	const long int a = 16807;
	const long int c = 0;
	const long m = 2147483647;
	static long int x0 = 17;
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

void Initialize (int nMol, double temperature, Vector* r, Vector* v)
{
	double x, y, z;
	DO_MOL
	{
		x = RandomGenerator(-0.5 * L, 0.5 * L);
		y = RandomGenerator(-0.5 * L, 0.5 * L);
		z = RandomGenerator(-0.5 * L, 0.5 * L);
		r[n].SetValues(x, y, z);
	}
	
	double sigmaMB = sqrt(temperature);
	double pMax = 1. / (sqrt(2.*M_PI) * temperature);
	double vRand, vTest, p;
	double a = -5. * sigmaMB;
	double b = 5. * sigmaMB;
	int molCount = 0;
	int coordCount = 0;
	double coordVal[3];
	
	while (molCount < nMol)
	{
		vRand = RandomGenerator(a, b);
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
}

void ComputeTotEnergy (int nMol, Vector* r, Vector* v, double & energy)
{
	double uSum = 0., vvSum = 0.;
	Vector dr;
	double ri6;
	DO_MOL vvSum += v[n] * v[n];
	for (int i=0; i<nMol-1; ++i)
	{
		for (int j=i+1; j<nMol; ++j)
		{
			dr = r[i] - r[j];
			PeriodicBC (dr);
			if (dr.Mag() < rMax)
			{
				ri6 = pow (1. / dr.Mag(),6);
				uSum += 4. * ri6 * (ri6 - 1.);
			}
		}
	}
	energy = (0.5 * vvSum + uSum) / nMol;
}

void ComputeStepEnergy (int nMol, int k, Vector rTest, Vector vTest, Vector* r, Vector* v, double & dE)
{
	Vector drN, drO;
	double localEOld = 0., localENew = 0.;
	double ri6N, ri6O;
	DO_MOL
	{
		if (n != k)
		{
			drN = rTest - r[n];
			PeriodicBC (drN);
			if (drN.Mag() < rMax)
			{
				drO = r[k] - r[n];
				PeriodicBC (drO);
				ri6O = pow (1. / drO.Mag(),6);
				ri6N = pow (1. / drN.Mag(),6);
				localEOld += 4. * ri6O * (ri6O - 1.);
				localENew += 4. * ri6N * (ri6N - 1.);
			}
		}
	}
	localEOld += 0.5 * (v[k] * v[k]);
	localENew += 0.5 * (vTest * vTest);
	dE = localENew - localEOld;	
}

void Metropolis (int nMol, double temperature, Vector* r, Vector* v)
{
	Vector dr, dv, rTest, vTest;
	double dE, ratio, limit;
	
	for (int k=0; k<nMol; ++k)
	{
		dr.SetValues (RAND_R, RAND_R, RAND_R);
		dv.SetValues (RAND_V, RAND_V, RAND_V);
		rTest = r[k] + dr;
		PeriodicBC (rTest);
		vTest = v[k] + dv;
		ComputeStepEnergy (nMol, k, rTest, vTest, r, v, dE);
		ratio = exp(-dE / temperature);
		if (ratio > 1.)
		{
			r[k] = rTest;
			v[k] = vTest;
			++ acceptCount;
		}
		else
		{
			limit = RandomGenerator(0., 1.);
			if (ratio >= limit)
			{
				r[k] = rTest;
				v[k] = vTest;
				++ acceptCount;
			}
		}
	++ totalCount;
	}	
}

void ComputePCF (int nMol, double rMaxPCF, double deltaR, int* histPCF, Vector* r)
{
	Vector dr;
	double R;
	int index;
	for (int i=0; i<nMol-1; ++i)
	{
		for (int j=i+1; j<nMol; ++j)
		{
			dr = r[i] - r[j];
			PeriodicBC (dr);
			R = dr.Mag();
			if (R < rMaxPCF)
			{
				index = (int) (R / deltaR);
				++ histPCF[index];
			}
		}
	}
}
void PlotResult (ostream & out, int stepCount, const double energy)
{
	out << stepCount << "\t" << energy << endl;
}

int main (int argc, char const* argv[])
{
	time_t tStart, tEnd;
	time (& tStart);
	const int stepLimit = 5000;
	const int startAccum = 100;
	double vol = L * L * L;
	
	const int histPCFSize = 100;
	double rMaxPCF = 0.5 * L;
	double deltaR = rMaxPCF / (double) histPCFSize;
	int* histPCF = new int[histPCFSize];
	for (int i=0; i<histPCFSize; ++i)
		histPCF[i] = 0;	
	
	double temperature, density;
	cout << "Temperatura: ";
	cin >> temperature;
	cout << "Densità: ";
	cin >> density;
	cout << "rRange: ";
	cin >> rRange;
	cout << "vRange: ";
	cin >> vRange;
	
	int nMol = (int) vol * density;
	Vector* r = new Vector [nMol];
	Vector* v = new Vector [nMol];
	Initialize (nMol, temperature, r, v);
	double energy = 0.;
	double eSum = 0., eeSum = 0.;
	double eMed, eErr;
	
	ofstream out1 ("Result.dat");
	out1.precision(10);
	ComputeTotEnergy (nMol, r, v, energy);
	int stepCount = 0;
	int propCount = 0;
	double acceptRatio;
	while (stepCount  < stepLimit)
	{
		++ stepCount;
		if (stepCount % 100 == 0) 
		{
			acceptRatio = (double) acceptCount / totalCount;
			cout << "step: " << stepCount  << " / " << stepLimit << "\tAcceptRatio: " << acceptRatio << endl;
		}
		
		Metropolis (nMol, temperature, r, v);
		if (stepCount > startAccum && stepCount % 100 == 0)
		{
			++ propCount;
			ComputeTotEnergy (nMol, r, v, energy);
			ComputePCF (nMol, rMaxPCF, deltaR, histPCF, r);
			eSum += energy;
			eeSum += energy * energy;
			PlotResult (out1, stepCount, energy);
		}
	}
	out1.close();
	delete [] r;
	delete [] v;
	
	ofstream out2 ("ResultPCF.dat");
	out2.precision (10);
	double histPCFDouble;
	for (int i=1; i<histPCFSize; ++i)
	{
		histPCFDouble = (double) histPCF[i] * 2. / ((double) (nMol - 1.) * density * pow(i * deltaR - 0.5*deltaR,2) * deltaR * 4.* M_PI * propCount);
		out2 << i*deltaR << "\t" << histPCFDouble << endl; 
	}
	out2.close();
	
	eMed = eSum / propCount;
	eErr = sqrt(((eeSum / propCount) - pow(eSum / propCount,2)) / propCount);
	cout << "In Result.dat: ";
	cout << "1: step | 2: EMed" << endl;
	cout << "Accept ratio: " << acceptRatio << endl;
	time (& tEnd);
	cout << "Energy: " << eMed << " +- " << eErr << endl;
	
	cout << "Done: execution time " << (tEnd - tStart) << " s\a" << endl;
	delete [] histPCF;	
	
	return 0;
}
