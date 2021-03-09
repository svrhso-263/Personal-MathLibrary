// MathLibrary.h
#pragma once

namespace MathLibrary {
	class Arithmatic {
	public:
		// Returns a + b.
		static double Add(double a, double b);

		// Returns a - b.
		static double Subtract(double a, double b);

		// Returns a * b.
		static double Multiply(double a, double b);

		// Returns a / b.
		static double Divide(double a, double b);

		// Returns a to the power of b.
		static float Power(float a, int b);

		// Returns the square root of a.
		static float SquareRoot(const float &a);

		// Returns true if a is an even number. 
		static bool IsEven(int a);

		// Returns true if a is an odd number.
		static bool IsOdd(int a);

		// Returns the absolute value of an number
		static float Absolute(float a);

		// Returns the midpoint of two numbers
		static float Midpoint(float a, float b);
	};

	class Equations {
	public:
		// Returns the length of the hypotenuse (opposite's length, adjacent's length). Returns as the length measurement you plugged in.
		static float HypoPythagoras(float opposite, float adjacent);

		// Returns the length of the adjacent or opposite (hypotenuse length, other sides length). Return as the length measurement you plugged in.
		static float OpAdPythagoras(float hypotenuse, float a);

		// Returns the attraction of two objects (gravitation force in m/s, mass1 in mass, mass2 in mass, r in distance ). Returns force in newtons.
		static float Gravity(float g, float mass1, float mass2, float r);

		// Returns the speed of an object falling towards the earth each second (mass in KG, surface area in m2, time in seconds). Returns gravity as m/s.
		static float EarthGravity(float mass, float surfaceArea, float time);

		// Returns the terminal velocity of an object (acceleration in m/s, mass in KG, surface area in m2). Returns terminal velocity in m/s.
		static float TerminalVelocity(float mass, float acceleration, float surfaceArea);

		// Returns the drag of an object (mass in KG, acceleration in m/s, time in seconds). Returns drag in newtons;
		static float Drag(float mass, float acceleration, float time);

		// Returns the degrees of an angle (float radians).
		static float RadianToDegrees(float rad);

		// Returns the radians of an angle (float degrees).
		static float DegreesToRadians(float deg);
	};

	class Trigonometry {
	private:
		// Rint used for calculating sine
		static float Rint(float a);

		// Cosine core used for calculating sine and cosine
		static float CosineCore(float angle);

		// Sine core used for calculating Sine
		static float SineCore(float angle);

		// Arc sine cor used for calculating sine
		static float ASineCore(float angle);

	public:
		// Returns the sine of an angle
		static float Sine(float angle);

		// Returns the cosine of an angle
		static float ArcCosine(float decimal);
	};
}