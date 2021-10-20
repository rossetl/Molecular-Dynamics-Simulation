#ifndef VECTOR2D_H_
#define VECTOR2D_H_
#include <cmath>

namespace VECTOR3D
{
	class Vector
	{
	private:
		double X;
		double Y;
		double Z;
	public:
		Vector();
		Vector(double n1, double n2, double n3);
		~Vector();
		inline double x () const {return X;};
		inline double y () const {return Y;};
		inline double z () const {return Z;};
		void SetValues (double n1, double n2, double n3);
		void SetX (double n);
		void SetY (double n);
		void SetZ (double n);
		double Mag ();
		double Vol ();
		Vector operator + (const Vector & b) const;
		Vector operator - (const Vector & b) const;
		Vector operator - () const;
		Vector operator * (double n) const;
		Vector operator / (double n) const;
		void operator += (const Vector & b);
		void operator -= (const Vector & b);
		friend Vector operator * (double n, const Vector & a);
		friend Vector operator / (double n, const Vector & a);
		friend double operator * (Vector & a, Vector & b);
	};
	
	//Methods
	Vector::Vector()
	{
		X = Y = Z = 0.;
	}
	
	Vector::Vector(double n1, double n2, double n3)
	{
		X = n1;
		Y = n2;
		Z = n3;
	}
	
	Vector::~Vector ()
	{
	}
	
	void Vector::SetValues (double n1, double n2, double n3)
	{
		X = n1;
		Y = n2;
		Z = n3;
	}
	
	void Vector::SetX (double n)
	{
		X = n;
	}
	
	void Vector::SetY (double n)
	{
		Y = n;
	}
	
	void Vector::SetZ (double n)
	{
		Z = n;
	}
	
	double Vector::Mag()
	{
		return sqrt(X*X + Y*Y + Z*Z);
	}
	
	double Vector::Vol ()
	{
		return X*Y*Z;
	}
	
	Vector Vector::operator + (const Vector & b) const
	{
		return Vector(X + b.X, Y + b.Y, Z + b.Z);
	}
	
	Vector Vector::operator - (const Vector & b) const
	{
		return Vector(X - b.X, Y - b.Y, Z - b.Z);
	}
	
	Vector Vector::operator - () const
	{
		return Vector(-X, -Y, -Z);
	}
	
	Vector Vector::operator * (double n) const
	{
		return Vector(n*X, n*Y, n*Z);
	}
	
	Vector Vector::operator / (double n) const
	{
		return Vector(X/n, Y/n, Z/n);
	}
	
	void Vector::operator += (const Vector & b)
	{
		X += b.X;
		Y += b.Y;
		Z += b.Z;
	}
	
	void Vector::operator -= (const Vector & b)
	{
		X -= b.X;
		Y -= b.Y;
		Z -= b.Z;
	}
	
	Vector operator * (double n, Vector & a)
	{
		return a*n;
	}
	
	double operator * (Vector & a, Vector & b)
	{
		return a.X*b.X + a.Y*b.Y + a.Z*b.Z;
	}













} //end namespace VECTOR3D
#endif
