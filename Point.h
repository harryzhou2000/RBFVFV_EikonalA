#ifndef _POINT_H
#define _POINT_H
#include "TypeDefine.h"
#include <cmath>
#include <iostream>

// defined data type: for 2d problems
namespace ScalarCfv
{
	class point2D
	{
	public:
		real x;
		real y;

		point2D(real xx, real yy)
			: x(xx), y(yy) {}

		point2D()
			: x(0), y(0) {}
		inline point2D operator/(real x) const;
		inline point2D &operator/=(real x);
		inline point2D operator*(real x) const;
		inline friend point2D operator*(real x, const point2D &);
		inline friend point2D operator+(const point2D &x, const point2D &y);
		inline friend point2D operator-(const point2D &x, const point2D &y);
		inline point2D &operator+=(const point2D &y);
		inline point2D &operator-=(const point2D &y);
		inline point2D &operator=(real x);
		inline real length() const { return sqrt(x * x + y * y); }

		inline real getOuterProduct(const point2D &v1, const point2D &v2);
		inline real getInnerProduct(const point2D &v1, const point2D &v2);
		// innerProduct, outerProduct
		inline point2D setZero(void) { return 0.0 * (*this); }
	};



	inline point2D point2D::operator/(real xx) const
	{
		return point2D(x / xx, y / xx);
	}
	inline point2D point2D::operator*(real xx) const
	{
		return point2D(x * xx, y * xx);
	}
	inline point2D operator*(real x, const point2D &r)
	{
		return point2D(x * r.x, x * r.y);
	}
	inline point2D &point2D::operator/=(real xx)
	{
		x /= xx;
		y /= xx;
		return *this;
	}
	inline point2D operator+(const point2D &x, const point2D &y)
	{
		return point2D(x.x + y.x, x.y + y.y);
	}
	inline point2D operator-(const point2D &x, const point2D &y)
	{
		return point2D(x.x - y.x, x.y - y.y);
	}
	inline point2D &point2D::operator+=(const point2D &yy)
	{
		x += yy.x;
		y += yy.y;
		return *this;
	}
	inline point2D &point2D::operator-=(const point2D &yy)
	{
		x -= yy.x;
		y -= yy.y;
		return *this;
	}

	inline point2D &point2D::operator=(real xx)
	{
		x = xx;
		y = xx;
		return *this;
	}

	inline real getOuterProduct(const point2D &v1, const point2D &v2)
	{
		return v1.x * v2.y - v2.x * v1.y;
	}

	inline real getInnerProduct(const point2D &v1, const point2D &v2)
	{
		return v1.x * v2.x + v1.y * v2.y;
	}

	typedef point2D point;
	typedef std::vector<point> pointVector;
	typedef std::vector<pointVector> pointVectorVector;
	std::ostream &operator<<(std::ostream &out, const point2D &p);
}


#endif