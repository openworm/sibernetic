/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

#ifndef _VectorMath_
#define _VectorMath_

#include <math.h>

const float PI=3.14159265f;


//An object to represent a 3D vector or a 3D point in space
class Vector3D
{
public:
	float x;									// the x value of this Vector3D
	float y;									// the y value of this Vector3D
	float z;									// the z value of this Vector3D

    Vector3D():  // Constructor to set x = y = z = 0
        x(0),
        y(0),
        z(0)
	{

	}

    Vector3D(float x, float y, float z):			// Constructor that initializes this Vector3D to the intended values of x, y and z
        x(x),
        y(y),
        z(z)
    {

	}

    Vector3D& operator= (const Vector3D &v)			// operator= sets values of v to this Vector3D. example: v1 = v2 means that values of v2 are set onto v1
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

    Vector3D operator+ (const Vector3D &v)	const			// operator+ is used to add two Vector3D's. operator+ returns a new Vector3D
	{
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

    Vector3D operator- (const Vector3D &v)	const			// operator- is used to take difference of two Vector3D's. operator- returns a new Vector3D
	{
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

    Vector3D operator* (float value)	const		// operator* is used to scale a Vector3D by a value. This value multiplies the Vector3D's x, y and z.
	{
		return Vector3D(x * value, y * value, z * value);
	}


    Vector3D operator/ (float value)	const		// operator/ is used to scale a Vector3D by a value. This value divides the Vector3D's x, y and z.
	{
		return Vector3D(x / value, y / value, z / value);
	}

    Vector3D& operator+= (const Vector3D &v)			// operator+= is used to add another Vector3D to this Vector3D.
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

    Vector3D& operator-= (const Vector3D &v)			// operator-= is used to subtract another Vector3D from this Vector3D.
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	Vector3D& operator*= (float value)			// operator*= is used to scale this Vector3D by a value.
	{
		x *= value;
		y *= value;
		z *= value;
		return *this;
	}

	Vector3D& operator/= (float value)			// operator/= is used to scale this Vector3D by a value.
	{
		x /= value;
		y /= value;
		z /= value;
		return *this;
	}


    Vector3D operator- ()	const					// operator- is used to set this Vector3D's x, y, and z to the negative of them.
	{
		return Vector3D(-x, -y, -z);
	}

    float length()		const						// length() returns the length of this Vector3D
	{
		return sqrt(x*x + y*y + z*z);
    }

	void unitize()								// unitize() normalizes this Vector3D that its direction remains the same but its length is 1.
	{
		float length = this->length();

		if (length == 0)
			return;

		x /= length;
		y /= length;
		z /= length;
	}

    Vector3D unit() const						// unit() returns a new Vector3D. The returned value is a unitized version of this Vector3D.
	{
		float length = this->length();

		if (length == 0)
			return *this;
		
		return Vector3D(x / length, y / length, z / length);
	}

    float scaleM(const Vector3D &vec) const
	{
		return x * vec.x + y * vec.y + z * vec.z;
	}

    float getLengthSq_fast() const
	{
		return x*x + y*y + z*z;
	}

    float operator*(const Vector3D &v) const
	{
		return x*v.x + y*v.y + z*v.z;
	}

    Vector3D operator%(const Vector3D &v) const
	{
		Vector3D t;
		t.x=y*v.z-z*v.y;
		t.y=z*v.x-x*v.z;
		t.z=x*v.y-y*v.x;
		return t;
	}

    bool operator==(const Vector3D &v) const
	{
		if((x==v.x)&&(y==v.y)&&(z==v.z))
			return true;
		else
			return false;
	}

    static Vector3D RotateVector1AroundVector2(const Vector3D &v1, const Vector3D &v2 ,float alpha)
	{
		if(v1==v2) return v1;

		Vector3D ort1,ort2,ort3;

		//alpha = 5;

		alpha*=(float)PI/180.0f;

		ort1 = v2.unit();
		ort2 = (v1%v2).unit();
		ort3 = (ort2%ort1).unit();
/*		
		Vector3D t1 = ort1*(v1*ort1);
		Vector3D t2 = ort2*(v1*ort3)*sin(alpha);
		Vector3D t3 = ort3*(v1*ort3)*cos(alpha);
		Vector3D r = t1+t2+t3;
*/
		return ort1*(v1*ort1) + ort2*(v1*ort3)*sin(alpha) + ort3*(v1*ort3)*cos(alpha);
	}
};






#endif
