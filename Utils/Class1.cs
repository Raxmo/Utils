using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Utils
{
	public static class Consts
	{
		public static double Circ = Math.Pow(Math.E, 2 * Math.PI);
	}

	public struct Vec
	{
		public double X;
		public double Y;

		public Vec(double x, double y)
		{
			X = x;
			Y = y;
		}

		public static Vec operator + (Vec a, Vec b)
		{
			return new Vec(a.X + b.X, a.Y + b.Y);
		}

		public static Vec operator - (Vec a, Vec b)
		{
			return new Vec(a.X - b.X, a.Y - b.Y);
		}

		public static Vec operator * (Vec a, double b)
		{
			return new Vec(a.X * b, a.Y * b);
		}

		public static Vec operator * (double a, Vec b)
		{
			return b * a;
		}

		public static Vec operator * (Vec a, Vec b) //complex multiplication
		{
			return new Vec(a.X * b.X - a.Y * b.Y, a.X * b.Y + a.Y * b.X);
		}

		public static Vec operator ^ (double a, Vec b)
		{
			return new Vec(Math.Cos(Math.Log(a) * b.Y), Math.Sin(Math.Log(a) * b.Y)) * Math.Pow(a, b.X);
		}
	}



}
