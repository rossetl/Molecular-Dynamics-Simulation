#ifndef VECTOR2D_H_
#define VECTOR2D_H_
#include <cmath>

namespace VECTOR2D
{
	class Vector
	{
	private:
		double X;
		double Y;
	public:
		Vector();
		Vector(double n1, double n2);
		~Vector();
		double x () const {return X;};
		double y () const {return Y;};
		void SetValues (double n1, double n2);
		void SetX (double n);
		void SetY (double n);
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
		X = Y = 0.;
	}
	
	Vector::Vector(double n1, double n2)
	{
		X = n1;
		Y = n2;
	}
	
	Vector::~Vector ()
	{
	}
	
	void Vector::SetValues (double n1, double n2)
	{
		X = n1;
		Y = n2;
	}
	
	void Vector::SetX (double n)
	{
		X = n;
	}
	
	void Vector::SetY (double n)
	{
		Y = n;
	}
	
	double Vector::Mag()
	{
		return sqrt(X*X + Y*Y);
	}
	
	double Vector::Vol ()
	{
		return X*Y;
	}
	
	Vector Vector::operator + (const Vector & b) const
	{
		return Vector(X + b.X, Y + b.Y);
	}
	
	Vector Vector::operator - (const Vector & b) const
	{
		return Vector(X - b.X, Y - b.Y);
	}
	
	Vector Vector::operator - () const
	{
		return Vector(-X, -Y);
	}
	
	Vector Vector::operator * (double n) const
	{
		return Vector(n*X, n*Y);
	}
	
	Vector Vector::operator / (double n) const
	{
		return Vector(X/n, Y/n);
	}
	
	void Vector::operator += (const Vector & b)
	{
		X += b.X;
		Y += b.Y;
	}
	
	void Vector::operator -= (const Vector & b)
	{
		X -= b.X;
		Y -= b.Y;
	}
	
	Vector operator * (double n, Vector & a)
	{
		return a*n;
	}
	
	double operator * (Vector & a, Vector & b)
	{
		return a.X*b.X + a.Y*b.Y;
	}













} //end namespace VECTOR
#endif
