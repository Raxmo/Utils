using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Drawing;

namespace Utils
{
	public static class Consts
	{
		public static double Circ = Math.Pow(Math.E, 2 * Math.PI);
	}

	public static class Functions
	{
		public static double Sigmoid(double x)
		{
			return 1 / (1 - Math.Exp(- x));
		}
	}

	public struct Vec
	{
		public double X;
		public double Y;

		public string Print()
		{
			return "X: " + X + ", Y:" + Y;
		}

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

		public static Vec operator / (Vec a, double b)
		{
			return new Vec(a.X / b, a.Y / b);
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

		public static double operator | (Vec a, Vec b)
		{
			return (a.X * b.X) + (a.Y * b.Y);
		}

		public static Vec operator ^ (Vec a, double b)
		{
			var pr = Math.Sqrt(a | a);
			var pt = Math.Atan2(a.Y, a.X);
			return (Math.Pow(pr, b) * (Math.E ^ (new Vec(0, b * pt))));
		}

		public static bool operator == (Vec a, Vec b)
		{
			return (a.X == b.X && a.Y == b.Y);
		}

		public static bool operator != (Vec a, Vec b)
		{
			return !(a == b);
		}

		public static explicit operator Vec(System.Drawing.Point p)
		{
			return new Vec(p.X, p.Y);
		}

		public static explicit operator Vec(System.Windows.Point p)
		{
			return new Vec(p.X, p.Y);
		}

	}

	public static class Matrix
	{
		public static double[,] SubM(double[,] mat, int row, int col)
		{
			var c = new double[mat.GetLength(0) - 1, mat.GetLength(1) - 1];
			
			for(int i = 0; i < c.GetLength(0); i++)
			{
				for(int j = 0; j < c.GetLength(1); j++)
				{
					int ioffset = 0;
					int joffset = 0;

					if(i >= row)
					{
						ioffset = 1;
					}
					if(j >= col)
					{
						joffset = 1;
					}

					c[i, j] = mat[i + ioffset, j + joffset];
				}
			}
			return c;
		}

		public static double Det(double[,] mat)
		{
			if(mat.GetLength(0) == 2 && mat.GetLength(1) == 2)
			{
				return mat[0, 0] * mat[1, 1] - mat[1, 0] * mat[0, 1];
			}
			else
			{
				double total = 0;
				for(int i = 0; i < mat.GetLength(0); i++)
				{
					total += Math.Pow(-1, i) * mat[i, 0] * Det(SubM(mat, i, 0));
				}
				return total;
			}
		}

		public static double[,] Add(double[,] a, double[,] b)
		{
			var c = new double[a.GetLength(0), a.GetLength(1)];
			for(int i = 0; i < a.GetLength(0); i++)
			{
				for (int j = 0; j < a.GetLength(1); j++)
				{
					c[i, j] = a[i, j] + b[i, j];
				}
			}
			return c;
		}

		public static double[,] Mult(double[,] a, double[,] b)
		{
			var c = new double[a.GetLength(0), b.GetLength(1)];

			for(int i = 0; i < c.GetLength(0); i++)
			{
				for(int j = 0; j < c.GetLength(1); j++)
				{
					double total = 0;

					for (int k = 0; k < b.GetLength(0); k++)
					{
						total += a[i, k] * b[k, j];
					}
					c[i, j] = total;
				}
			}

			return c;
		}

		public static double[,] Mult(double[,] a, double b)
		{
			for (int i = 0; i < a.GetLength(0); i++)
			{
				for (int j = 0; j < a.GetLength(1); j++)
				{
					a[i, j] *= b;
				}
			}
			return a;
		}

		public static double[,] Mult(double a, double[,] b)
		{
			return Mult(b, a);
		}
	}

}
