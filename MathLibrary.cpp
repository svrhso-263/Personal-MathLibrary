// MathLibrary.cpp
// compile with: cl /c /EHsc MathLibrary.cpp
// post-build command: lib MathLibrary.obj

#include <cmath>
#include "MathLibrary.h"

namespace MathLibrary {

	double Arithmatic::Add(double a, double b) {
		return a + b;
	}

	double Arithmatic::Subtract(double a, double b) {
		return a - b;
	}

	double Arithmatic::Multiply(double a, double b) {
		return a * b;
	}

	double Arithmatic::Divide(double a, double b) {
		return a / b;
	}

	float Arithmatic::Power(float a, int b)
	{
		int i, result = 1;

		for (i = 0; i < b; i++) {
			result = result * a;
		}

		return result;
	}

	float Arithmatic::SquareRoot(const float &a)
	{
		static union { int i; float f; } u;
		u.i = 0x5F375A86 - (*(int*)&a >> 1);
		return (int(3) - a * u.f * u.f) * a * u.f * 0.5f;
	}

	bool Arithmatic::IsEven(int a)
	{
		if ((a % 2) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	bool Arithmatic::IsOdd(int a) {
		if ((a % 2) != 0) {
			return true;
		}
		else {
			return false;
		}
	}

	float Arithmatic::Absolute(float a) {
		return a * ((a > 0) - (a < 0));
	}

	float Arithmatic::Midpoint(float a, float b)
	{
		float ab = a + b;

		return ab / 2;
	}

	float Equations::HypoPythagoras(float opposite, float adjacent)
	{
		const float &squareHypo = (opposite * opposite) + (adjacent * adjacent);

		static union { int i; float f; } u;
		u.i = 0x5F375A86 - (*(int*)&squareHypo >> 1);
		float hypotenuse = (int(3) - squareHypo * u.f * u.f) * squareHypo * u.f * 0.5f;

		return hypotenuse;
	}

	float Equations::OpAdPythagoras(float hypotenuse, float a)
	{
		const float& squareB = (hypotenuse * hypotenuse) - (a * a);

		static union { int i; float f; } u;
		u.i = 0x5F375A86 - (*(int*)&squareB >> 1);
		float b = (int(3) - squareB * u.f * u.f) * squareB * u.f * 0.5f;

		return b;
	}

	float Equations::Gravity(float g, float m1, float m2, float r)
	{
		float force = g *((m1 * m2) / (r * r));

		return force;
	}

	float Equations::EarthGravity(float mass, float surfaceArea, float time) {
		float velocity = (9.8 * time) / 2;

		float terminalVelocity = MathLibrary::Arithmatic::SquareRoot((2 * mass * 9.8) / (1.225 * surfaceArea * 0.8));
		if (velocity > terminalVelocity)
			velocity = terminalVelocity;

		return velocity;
	}

	float Equations::TerminalVelocity(float mass, float acceleration, float surfaceArea) {
		float terminalVelocity = MathLibrary::Arithmatic::SquareRoot((2 * mass * acceleration) / (1.225 * surfaceArea * 0.8));

		return terminalVelocity;
	}

	float Equations::Drag(float mass, float acceleration, float time) {
		float velocity = acceleration * time;
		float drag = 0.8 * ((mass * MathLibrary::Arithmatic::Power(velocity, 2)) / 2) * acceleration;

		return drag;
	}

	float Equations::RadianToDegrees(float rad)
	{
		float deg = rad * (180 / 3.14159265359);

		return deg;
	}

	float Equations::DegreesToRadians(float deg)
	{
		float rad = deg * (3.14159265359 / 180);

		return rad;
	}

	float Trigonometry::Rint(float a)
	{
		float t = floor(fabs(a) + 0.5);
		return (a < 0.0) ? -t : t;
	}
	float Trigonometry::CosineCore(float angle)
	{
		float angle8, angle4, angle2;
		angle2 = angle * angle;
		angle4 = angle2 * angle2;
		angle8 = angle4 * angle4;
		return (-2.7236370439787708e-7 * angle2 + 2.4799852696610628e-5) * angle8 +
			(-1.3888885054799695e-3 * angle2 + 4.1666666636943683e-2) * angle4 +
			(-4.9999999999963024e-1 * angle2 + 1.0000000000000000e+0);
	}
	float Trigonometry::SineCore(float angle)
	{
		float angle4, angle2, t;
		angle2 = angle * angle;
		angle4 = angle2 * angle2;
		return ((2.7181216275479732e-6 * angle2 - 1.9839312269456257e-4) * angle4 +
			(8.3333293048425631e-3 * angle2 - 1.6666666640797048e-1)) * angle2 * angle + angle;
	}
	float Trigonometry::ASineCore(float angle)
	{
		float angle8, angle4, angle2;
		angle2 = angle * angle;
		angle4 = angle2 * angle2;
		angle8 = angle4 * angle4;
		return (((4.5334220547132049e-2 * angle2 - 1.1226216762576600e-2) * angle4 +
			(2.6334281471361822e-2 * angle2 + 2.0596336163223834e-2)) * angle8 +
			(3.0582043602875735e-2 * angle2 + 4.4630538556294605e-2) * angle4 +
			(7.5000364034134126e-2 * angle2 + 1.6666666300567365e-1)) * angle2 * angle + angle;
	}
	float Trigonometry::Sine(float angle) {
		float q, t;
		int quadrant;
		q = Rint(angle * 6.3661977236758138e-1);
		quadrant = (int)q;
		t = angle - q * 1.5707963267923333e+00;
		t = t - q * 2.5633441515945189e-12;
		if (quadrant & 1) {
			t = CosineCore(t);
		}
		else {
			t = SineCore(t);
		}
		return (quadrant & 2) ? -t : t;
	}
	float Trigonometry::ArcCosine(float decimal) {
		float xa, t;
		xa = fabs(decimal);

		if (xa > 0.5625) {
			t = 2.0 * ASineCore(MathLibrary::Arithmatic::SquareRoot(0.5 * (1.0 - xa)));
		}
		else {
			t = 1.5707963267948966 - ASineCore(xa);
		}

		return (decimal < 0.0) ? (3.1415926535897932 - t) : t;
	}
}