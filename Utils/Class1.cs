using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Drawing;
using System.Collections;

namespace Utils
{
	/* mathie utils framework and notes
	 * 
	 * scalars :			x
	 *	imagionary unit :	i
	 *	
	 * vectors :			{x, y}
	 *	complex numbers :	{r, i}
	 *	polars :			{r, Θ}
	 *	3polars :			{r, h, a}
	 *	quaterneons :		{r, i, j, k}
	 *	
	 * matricies:			{{00, 01},
	 *						 {10, 11}}
	 *			  
	 * 
	 * Tensors :
	 *		All tensors are arrays of tensors
	 * 
	 */

	public static class Consts
	{
		public static readonly double Circ = Math.Exp(2 * Math.PI);

		public static readonly Complex i = new Complex(0, 1);

		public static Random rand = new Random();

		public static double nrand(double mean, double sd)
		{
			var u1 = 1.0 - rand.NextDouble();
			var u2 = 1.0 - rand.NextDouble();

			var randNorm = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2 * Math.PI * u2);

			return randNorm * sd + mean;
		}
		


		public static double Sigmoid(double x)
		{
			return 1 / (1 + Math.Exp(-x));
		}
		
		public static Vec Sigmoid(Vec x)
		{
			Vec v = new Vec();
			foreach (double d in v)
			{
				v.Add(Sigmoid(d));
			}
			return v;
		}

		public static Matrix Sigmoid(Matrix x)
		{
			Matrix m = new Matrix();
			foreach (Vec v in x)
			{
				m.Add(Sigmoid(v));
			}
			return m;
		}



		public static double SigmoidN(double x, double mean, double sd, double scale)
		{
			return Sigmoid((x - mean) / sd) * scale;
		}

		public static Vec SigmoidN(Vec x, double mean, double sd, double scale)
		{
			Vec v = new Vec();
			foreach (double d in x)
			{
				v.Add(SigmoidN(d, mean, sd, scale));
			}
			return v;
		}

		public static Matrix SigmoidN(Matrix x, double mean, double sd, double scale)
		{
			Matrix m = new Matrix();
			foreach(Vec v in x)
			{
				m.Add(SigmoidN(v, mean, sd, scale));
			}
			return m;
		}



		public static Vec SoftMax(Vec x)
		{
			Vec s = new Vec();
			double sum = 0;
			foreach(double d in x)
			{
				sum += Math.Exp(d);
			}
			foreach(double d in x)
			{
				s.Add(Math.Exp(d) / sum);
			}
			return s;
		}
	}

	// Tensor stuffs

	public class Tensor : IEnumerable
	{
		// Tensor variable

		private List<dynamic> t = new List<dynamic>();
		private int depth = 0;

		// Tensor enumerator/indexing stuffs

		public IEnumerator GetEnumerator()
		{
			yield return t;
		}

		public void Add(params double[] a)
		{
			if (a.Length == 1)
			{
				t.Add(a[0]);
				depth = 1;
			}
			else
			{
				Tensor e = new Tensor();

				for (int el = 0; el < a.Length; el++)
				{
					e.Add(a[el]);
				}

				t.Add(e);
				depth = 2;
			}
		}

		public dynamic this[int index]
		{
			get
			{
				return t[index];
			}
			set
			{
				t[index] = value;
			}
		}
		
		// Tensor constructors
				
		// Tensor utilities

		public int Depth
		{
			get => depth;
		}

		public int Level
		{
			get
			{
				return t.Count();
			}
		}

		public Tensor Copy()
		{
			var c = new Tensor();
			c.depth = depth;
			for (int spot = 0; spot < t.Count; spot++)
			{
				if (t[spot].GetType() == typeof(Tensor))
				{
					c.t.Add(t[spot].Copy());
				}
				else
				{
					c.t.Add(t[spot]);
				}
			}
			return c;
		}

		// Tensor operations

		public static Tensor operator + (Tensor a, Tensor b)
		{
			var c = new Tensor();
			for (int i = 0; i < a.Level; i++)
			{
				c.t.Add(a[i] + b[i]);
			}
			return c;
		}

		public static Tensor operator - (Tensor  a, Tensor b)
		{
			var c = new Tensor();
			for (int i = 0; i < a.Level; i++)
			{
				c.t.Add(a[i] - b[i]);
			}
			return c;
		}

		public static Tensor operator * (Tensor a, double b)
		{
			var c = new Tensor();
			for (int i = 0; i < a.Level; i++)
			{
				c.t.Add(a[i] * b);
			}
			return c;
		}

		public static Tensor operator * (double a, Tensor b)
		{
			return b * a;
		}

		public static Tensor operator / (Tensor a, double b)
		{
			var c = new Tensor();
			for (int i = 0; i < a.Level; i++)
			{
				c.t.Add(a[i] / b);
			}
			return c;
		}
	}

	// Vector stuffs

	public class Vec : IEnumerable
	{
		// vector variables

		private List<double> v = new List<double>();

		// Vector enumeration and indexing

		public IEnumerator GetEnumerator()
		{
			yield return v;
		}

		public void Add(params double[] a)
		{
			v.Add(a[0]);
		}

		public void Add(double d)
		{
			v.Add(d);
		}

		public double this[int dim]
		{
			get
			{
				return v[dim];
			}
			set
			{
				v[dim] = value;
			}
		}

		// Vector constructors

		public Vec(int dims)
		{
			for (int dim = 0; dim < dims; dim++)
			{
				v.Add(0);
			}
		}

		public Vec()
		{

		}

		public Vec(double[] a)
		{
			v = a.ToList();
		}

		public static explicit operator Vec(double[] a)
		{
			return new Vec(a);
		}

		// Vector utilities

		public static Vec operator ~ (Vec a)
		{
			return a / a.Mag;
		}

		public string Print()
		{
			string p = "[";
			foreach (double d in v)
			{
				p += " " + d;
			}
			p += " ]";
			return p;
		}
		
		public int Dims
		{
			get
			{
				return v.Count();
			}
		}

		public Vec Copy()
		{
			var c = new Vec();
			foreach (double d in v)
			{
				c.v.Add(d);
			}
			return c;
		}

		public double Mag
		{
			get => Math.Sqrt(this * this);
		}

		public void RemoveDim(int dim)
		{
			v.RemoveAt(dim);
		}

		public void RemoveDims(int start, int num)
		{
			v.RemoveRange(start, num);
		}

		public void Concat(Vec a)
		{
			v.AddRange(a.v);
		}

		// Vector operations

		public static Vec operator + (Vec a, Vec b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c.v.Add(a[dim] + b[dim]);
			}
			return c;
		}

		public static Vec operator - (Vec a, Vec b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c.v.Add(a[dim] - b[dim]);
			}
			return c;
		}

		public static Vec operator * (Vec a, double b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c.v.Add(a[dim] * b);
			}
			return c;
		}

		public static Vec operator * (double a, Vec b)
		{
			return b * a;
		}

		public static Vec operator / (Vec a, double b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c.v.Add(a[dim] / b);
			}
			return c;
		}

		public static double operator * (Vec a, Vec b)
		{
			double sum = 0;
			for (int dim = 0; dim < a.Dims; dim++)
			{
				sum += a[dim] * b[dim];
			}
			return sum;
		}

		// Vector convertions

		public static explicit operator Complex (Vec a)
		{
			return new Complex(a.v.ToArray());
		}

		public static explicit operator Quaterneon (Vec a)
		{
			return new Quaterneon(a.v.ToArray());
		}

		public static explicit operator Polar(Vec a)
		{
			Polar p = new Polar();
			p.r = a.Mag;
			p.t = Math.Atan2(a[0], a[1]);
			return p;
		}

		public static explicit operator Polar3(Vec a)
		{
			Polar3 p = new Polar3();
			p.r = a.Mag;
			p.h = Math.Acos(a[2] / p.r);
			p.a = Math.Atan2(a[1], a[0]);
			return p;
		}
	}

	// Polar stuffs

	public class Polar
	{
		// polar variables

		private double[] p = new double[2];

		// polar constructors

		public Polar() { }

		public Polar(double r, double Θ)
		{
			this.r = r;
			t = Θ;
		}

		// polar getters and such

		public double t
		{
			get => p[1];
			set => p[1] = value;
		}

		public double r
		{
			get => p[0];
			set => p[0] = value;
		}

		// polar utils

		public Polar Copy()
		{
			return new Polar(r, t);
		}

		public static Polar operator ~ (Polar a)
		{
			return new Polar(1, a.t);
		}


		// polar operations

		public static Polar operator + (Polar a, Polar b)
		{
			return new Polar(a.r + b.r, a.t + b.t);
		}

		public static Polar operator - (Polar a, Polar b)
		{
			return new Polar(a.r - b.r, a.t - b.t);
		}

		public static Polar operator * (Polar a, double b)
		{
			return new Polar(a.r * b, a.t);
		}

		public static Polar operator * (double b, Polar a)
		{
			return a * b;
		}

		public static Polar operator / (Polar a, double b)
		{
			return new Polar(a.r / b, a.t);
		}

		// polar conversions
		
		public static explicit operator Vec(Polar a)
		{
			Vec v = new Vec();
			v.Add(a.r * Math.Cos(a.t));
			v.Add(a.r * Math.Sin(a.t));
			return v;
		}
	}

	// polar 3 stuffs

	public class Polar3
	{
		// polar 3 variables

		private double[] p3 = new double[3];

		// polar 3 constructors

		public Polar3() { }

		public Polar3(double r, double h, double a)
		{
			p3 = new double[3] { r, h, a };
		}

		// polar 3 getters/setters

		public double r
		{
			get => p3[0];
			set => p3[0] = value;
		}

		public double h
		{
			get => p3[1];
			set => p3[1] = value;
		}

		public double a
		{
			get => p3[2];
			set => p3[2] = value;
		}

		// polar 3 utils

		// polar 3 operations

		public static Polar3 operator + (Polar3 a, Polar3 b)
		{
			return new Polar3(a.r + b.r, a.h + b.h, a.a + b.a);
		}

		public static Polar3 operator - (Polar3 a, Polar3 b)
		{
			return new Polar3(a.r - b.r, a.h - b.h, a.a - b.a);
		}

		public static Polar3 operator * (Polar3 a, double b)
		{
			return new Polar3(a.r * b, a.h, a.a);
		}

		public static Polar3 operator * (double b, Polar3 a)
		{
			return a * b;
		}

		public static Polar3 operator / (Polar3 a, double b)
		{
			return new Polar3(a.r / b, a.h, a.a);
		}

		// polar 3 conversions

		public static explicit operator Vec(Polar3 a)
		{
			Vec v = new Vec();
			v.Add(a.r * Math.Sin(a.h) * Math.Cos(a.a));
			v.Add(a.r * Math.Sin(a.h) * Math.Sin(a.a));
			v.Add(a.r * Math.Cos(a.h));
			return v;
		}

		public static explicit operator Polar(Polar3 a)
		{
			return new Polar(a.a, a.h);
		}
	}

	// Matrix stuffs

	public class Matrix : IEnumerable
	{
		// Matrix Variables

		private List<Vec> m = new List<Vec>();

		// Matrix enumeration and indexing

		public IEnumerator GetEnumerator()
		{
			yield return m;
		}

		public void Add(params double[] a)
		{
			m.Add(new Vec(a));
		}

		public void Add(Vec v)
		{
			m.Add(v);
		}
				
		public dynamic this[int row, int col]
		{
			get
			{
				if (row < 0 && col >= 0)
				{
					Vec c = new Vec();
					for (int r = 0; r < Rows; r++)
					{
						c.Add(m[r][col]);
					}
					return c;
				}
				if (col < 0 && row >= 0)
				{
					return m[row];
				}
				if (row >= 0 && col >= 0)
				{
					return m[row][col];
				}
				else
				{
					return null;
				}
			}
			set
			{
				if (col < 0)
				{
					m[row] = value;
				}
				if (row < 0)
				{
					for (int r = 0; r < Rows; r++)
					{
						m[r][col] = value[r];
					}
				}
				if (row >=0 && col >= 0)
				{
					m[row][col] = value;
				}
			}
		}

		// matrix constructors

		public Matrix() { }
		
		public Matrix(int rows, int cols)
		{
			for (int row = 0; row < rows; row++)
			{
				m.Add(new Vec(cols));
			}
		}

		// matrix utilities

		public string Print()
		{
			string p = "\n[";
			for (int row = 0; row < Rows; row++)
			{
				if (row > 0)
				{
					p += " ";
				}
				p += m[row].Print();
				if (row == Rows - 1)
				{
					p += "]\n";
				}
				else
				{
					p += "\n";
				}
			}
			return p;
		}

		public Matrix Copy()
		{
			Matrix c = new Matrix();
			for (int row = 0; row < Rows; row++)
			{
				c.m.Add(m[row].Copy());
			}
			return c;
		}

		public int Rows
		{
			get
			{
				return m.Count();
			}
		}

		public int Cols
		{
			get
			{
				return m[0].Dims;
			}
		}

		public double Det
		{
			get
			{
				if (Rows == Cols)
				{
					if (Rows == 2)
					{
						return this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];
					}
					else
					{
						double sum = 0;
						for (int col = 0; col < Cols; col++)
						{
							double fac = Math.Pow(-1, col % 2);
							sum += this[0, col] * fac * MinorMatrix(0, col).Det;
						}
						return sum;
					}
				}
				else
				{
					return double.NaN;
				}
			}
		}

		public static Matrix Identity(int size)
		{
			var i = new Matrix(size, size);
			for (int place = 0; place < size; place++)
			{
				i[place, place] = 1;
			}
			return i;
		}

		public static Matrix Rotation2D(double Θ)
		{
			return new Matrix() { { Math.Cos(Θ), -Math.Sin(Θ) }, { Math.Sin(Θ), Math.Cos(Θ) } };
		}

		public Matrix Augment(Matrix b)
		{
			Matrix aug = Copy();
			for (int row = 0; row < Rows; row++)
			{
				aug.m[row].Concat(b.m[row]);
			}
			return aug;
		}

		public void DeAugment(int cols)
		{
			for (int row = 0; row < Rows; row++)
			{
				m[row].RemoveDims(0, cols);
			}
		}

		public Matrix MinorMatrix(int row, int col)
		{
			Matrix mm = Copy();
			if (row >= 0)
			{
				mm.m.RemoveAt(row);
			}
			if (col >= 0)
			{
				for (int r = 0; r < mm.Rows; r++)
				{
					mm.m[r].RemoveDim(col);
				}
			}
			return mm;
		}

		public Matrix Inverse()
		{
			Matrix inv = Copy();
			inv = inv.Augment(Identity(inv.Rows));

			// setting things to identity

			for (int col = 0; col < Rows; col++)
			{
				inv[col, -1] /= inv[col, col];

				for (int row = 0; row < Rows; row++)
				{
					if (row != col)
					{
						inv[row, -1] -= inv[col, -1] * inv[row, col];
					}
				}
			}
			inv.DeAugment(Rows);
			return inv;
		}
		
		public Matrix Transpose()
		{
			var t = new Matrix();
			for (int col = 0; col < Cols; col++)
			{
				t.m.Add(this[-1, col]);
			}
			return t;
		}

		// matrix operations

		public static Matrix operator + (Matrix a, Matrix b)
		{
			Matrix c = new Matrix();
			for (int row = 0; row < a.Rows; row++)
			{
				c.m.Add(a[row, -1] + b[row, -1]);
			}
			return c;
		}

		public static Matrix operator - (Matrix a, Matrix b)
		{
			Matrix c = new Matrix();
			for (int row = 0; row < a.Rows; row++)
			{
				c.m.Add(a[row, -1] - b[row, -1]);
			}
			return c;
		}

		public static Matrix operator * (Matrix a, double b)
		{
			Matrix s = new Matrix();
			for (int row = 0; row < a.Rows; row++)
			{
				s.Add(a[row, -1] * b);
			}
			return s;
		}

		public static Matrix operator * (double b, Matrix a)
		{
			return a * b;
		}

		public static Matrix operator * (Matrix a, Matrix b)
		{
			if (a.Cols == b.Rows)
			{
				Matrix g = new Matrix();
				for (int row = 0; row < a.Rows; row++)
				{
					Vec r = new Vec();
					for (int col = 0; col < b.Cols; col++)
					{
						r.Add(a[row, -1] * b[-1, col]);
					}
					g.Add(r);
				}
				return g;
			}
			else
			{
				return null;
			}
		}

		public static Matrix operator / (Matrix a, double b)
		{
			Matrix s = new Matrix();
			for (int row = 0; row < a.Rows; row++)
			{
				s.Add(a[row, -1] / b);
			}
			return s;
		}

		public static Vec operator * (Matrix a, Vec b)
		{
			if (a.Cols == b.Dims)
			{
				Vec p = new Vec();
				for (int r = 0; r < b.Dims; r++)
				{
					p.Add(a[r, -1] * b);
				}
				return p;
			}
			else
			{
				return null;
			}
		}
	}

	// Complex number stuffs

	public class Complex
	{
		// Complex number variables

		private double[] z = new double[2];

		// Complex getters and setters

		public double r
		{
			get
			{
				return z[0];
			}
			set
			{
				z[0] = value;
			}
		}

		public double i
		{
			get
			{
				return z[1];
			}
			set
			{
				z[1] = value;
			}
		}

		public double Abs
		{
			get => new Vec(z).Mag;
		}

		// Complex constructors

		public Complex() { }

		public Complex(double[] c)
		{
			z = c;
		}

		public Complex(double r, double i)
		{
			this.r = r;
			this.i = i;
		}

		// Complex operations

		public static Complex operator + (Complex a, Complex b)
		{
			Complex c = new Complex();
			c.r = a.r + b.r;
			c.i = a.i + b.i;
			return c;
		}

		public static Complex operator - (Complex a, Complex b)
		{
			Complex c = new Complex();
			c.r = a.r - b.r;
			c.i = a.i - b.i;
			return c;
		}

		public static Complex operator * (Complex a, double b)
		{
			Complex c = new Complex();
			c.r = a.r * b;
			c.i = a.i * b;
			return c;
		}

		public static Complex operator * (double b, Complex a)
		{
			return a * b;
		}

		public static Complex operator * (Complex a, Complex b)
		{
			Complex c = new Complex();
			c.r = a.r * b.r - a.i * b.i;
			c.i = a.r * b.i + a.i * b.r;
			return c;
		}

		public static Complex operator / (Complex a, double b)
		{
			Complex c = new Complex();
			c.r = a.r / b;
			c.i = a.i / b;
			return c;
		}

		public static Complex operator ^ (double a, Complex b)
		{
			return new Complex(Math.Cos(Math.Log(a) * b.i), Math.Sin(Math.Log(a) * b.i)) * Math.Pow(a, b.r);
		}

		// Complex Functions

		// complex convertions

		public static explicit operator Vec(Complex a)
		{
			return new Vec(a.z);
		}
	}

	// Quaterneon stuffs

	public class Quaterneon
	{
		// Quaterneon variables

		private double[] q = new double[4];

		// Quaterneon getters and such

		public double r
		{
			get
			{
				return q[0];
			}
			set
			{
				q[0] = value;
			}
		}

		public double i
		{
			get
			{
				return q[1];
			}
			set
			{
				q[1] = value;
			}
		}

		public double j
		{
			get => q[2];
			set => q[2] = value;
		}

		public double k
		{
			get => q[3];
			set => q[3] = value;
		}

		public double Abs
		{
			get => new Vec(q).Mag;
		}

		// Quaterneon constructors

		public Quaterneon() { }
		
		public Quaterneon(double r, double i, double j, double k)
		{
			this.r = r;
			this.i = i;
			this.j = j;
			this.k = k;
		}

		public Quaterneon(double[] q)
		{
			this.q = q;
		}

		public Quaterneon(double s, Vec v)
		{
			r = s;
			i = v[0];
			j = v[1];
			k = v[2];
		}

		// Quaterneon utils

		public string Print()
		{
			return ("< " + r + " " + i + " " + j + " " + k + " >");
		}

		public static explicit operator Vec(Quaterneon q)
		{
			return new Vec(q.q);
		}

		public Vec ToV3()
		{
			return new Vec(new double[3] { i, j, k });
		}

		/// <summary>
		/// Creates a rotation quaterneon from an angle and an axis vector
		/// </summary>
		/// <param name="Θ">Theta of rotation</param>
		/// <param name="v">Axis of rotation</param>
		/// <returns>Rotation quaterneon</returns>
		public static Quaterneon Rotation(double Θ, Vec v)
		{
			return new Quaterneon(Math.Cos(Θ / 2), (~v) * Math.Sin(Θ / 2));
		}

		// Quaterneon operations

		public static Quaterneon operator + (Quaterneon a, Quaterneon b)
		{
			return (Quaterneon)((Vec)a + (Vec)b);
		}

		public static Quaterneon operator - (Quaterneon a, Quaterneon b)
		{
			return (Quaterneon)((Vec)a - (Vec)b);
		}

		public static Quaterneon operator * (Quaterneon a, double b)
		{
			return (Quaterneon)((Vec)a * b);
		}

		public static Quaterneon operator * (double b, Quaterneon a)
		{
			return a * b;
		}

		public static Quaterneon operator / (Quaterneon a, double b)
		{
			return (Quaterneon)((Vec)a / b);
		}

		public static Quaterneon operator * (Quaterneon a, Quaterneon b)
		{
			Quaterneon c = new Quaterneon();
			c.r = a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k;
			c.i = a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j;
			c.j = a.r * b.j - a.i * b.k + a.j * b.r + a.k * b.i;
			c.k = a.r * b.k + a.i * b.j + a.j * b.i + a.k * b.r;
			return c;
		}

		/// <summary>
		/// Congugation of a Quaterneon
		/// </summary>
		/// <param name="a">Base Quaterneon</param>
		/// <returns>a new quaterneon that is a congugate of the root quaterneon</returns>
		public static Quaterneon operator ! (Quaterneon a)
		{
			return new Quaterneon(a.r, -a.i, -a.j, -a.k);

		}

		/// <summary>
		/// Rotational multiplication of quaterneons
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <returns></returns>
		public static Quaterneon operator ^ (Quaterneon a, Quaterneon b)
		{
			return b * a * (!b);
		}
		
		// Quaterneon functions
	}
	
	//public static class Consts
	//{
	//	public static double Circ = Math.Pow(Math.E, 2 * Math.PI);
	//}

	//public static class Functions
	//{
	//	public static double Sigmoid(double x)
	//	{
	//		return 1 / (1 - Math.Exp(- x));
	//	}
	//}

	//public struct Vec
	//{
	//	public double X;
	//	public double Y;

	//	public string Print()
	//	{
	//		return "X: " + X + ", Y:" + Y;
	//	}

	//	public Vec(double x, double y)
	//	{
	//		X = x;
	//		Y = y;
	//	}

	//	public static Vec operator + (Vec a, Vec b)
	//	{
	//		return new Vec(a.X + b.X, a.Y + b.Y);
	//	}

	//	public static Vec operator - (Vec a, Vec b)
	//	{
	//		return new Vec(a.X - b.X, a.Y - b.Y);
	//	}

	//	public static Vec operator * (Vec a, double b)
	//	{
	//		return new Vec(a.X * b, a.Y * b);
	//	}

	//	public static Vec operator / (Vec a, double b)
	//	{
	//		return new Vec(a.X / b, a.Y / b);
	//	}

	//	public static Vec operator * (double a, Vec b)
	//	{
	//		return b * a;
	//	}

	//	public static Vec operator * (Vec a, Vec b) //complex multiplication
	//	{
	//		return new Vec(a.X * b.X - a.Y * b.Y, a.X * b.Y + a.Y * b.X);
	//	}

	//	public static Vec operator ^ (double a, Vec b)
	//	{
	//		return new Vec(Math.Cos(Math.Log(a) * b.Y), Math.Sin(Math.Log(a) * b.Y)) * Math.Pow(a, b.X);
	//	}

	//	public static double operator | (Vec a, Vec b)
	//	{
	//		return (a.X * b.X) + (a.Y * b.Y);
	//	}

	//	public static Vec operator ^ (Vec a, double b)
	//	{
	//		var pr = Math.Sqrt(a | a);
	//		var pt = Math.Atan2(a.Y, a.X);
	//		return (Math.Pow(pr, b) * (Math.E ^ (new Vec(0, b * pt))));
	//	}

	//	public static bool operator == (Vec a, Vec b)
	//	{
	//		return (a.X == b.X && a.Y == b.Y);
	//	}

	//	public static bool operator != (Vec a, Vec b)
	//	{
	//		return !(a == b);
	//	}

	//	public static explicit operator Vec(System.Drawing.Point p)
	//	{
	//		return new Vec(p.X, p.Y);
	//	}

	//	public static explicit operator Vec(System.Windows.Point p)
	//	{
	//		return new Vec(p.X, p.Y);
	//	}

	//	public static Vec operator %(Vec a, Vec b)
	//	{
	//		return (a + b) / 2;
	//	}

	//}

	//public class Matrix : IEnumerable
	//{
	//	/* What we know abuot matricies
	//	 * - matricies are defined by rows, cols
	//	 * - addition is element-wise
	//	 */

	//	public IEnumerator GetEnumerator()
	//	{
	//		yield return m;
	//	}

	//	private int currentRow = 0;
	//	public void Add(params double[] a)
	//	{
	//		for (int c = 0; c < a.Length; c++)
	//		{
	//			m[currentRow, c] = a[c];
	//		}
	//		currentRow++;
	//	}

	//	private double[,] m;

	//	public double this[int rows, int cols]
	//	{
	//		get => m[rows, cols];
	//		set => m[rows, cols] = value;
	//	}		

	//	public Matrix(int rows, int cols)
	//	{
	//		m = new double[rows, cols];
	//	}

	//	public Matrix(Matrix a)
	//	{
	//		m = a.m;
	//	}

	//	public int Rows()
	//	{
	//		return m.GetLength(0);
	//	}

	//	public int Cols()
	//	{
	//		return m.GetLength(1);
	//	}

	//	public static Matrix operator + (Matrix a, Matrix b)
	//	{
	//		var c = new Matrix(a.Rows(), a.Cols());

	//		for (int rows = 0; rows < a.Rows(); rows++)
	//		{
	//			for (int cols = 0; cols < a.Cols(); cols++)
	//			{
	//				c[rows, cols] = a[rows, cols] + b[rows, cols];
	//			}
	//		}

	//		return c;
	//	}

	//	public static Matrix operator - (Matrix a, Matrix b)
	//	{
	//		var c = new Matrix(a.Rows(), a.Cols());

	//		for (int rows = 0; rows < a.Rows(); rows++)
	//		{
	//			for (int cols = 0; cols < a.Cols(); cols++)
	//			{
	//				c[rows, cols] = a[rows, cols] - b[rows, cols];
	//			}
	//		}

	//		return c;
	//	}

	//	public static Matrix operator * (Matrix a, double b)
	//	{
	//		var c = a;

	//		for (int rows = 0; rows < a.Rows(); rows++)
	//		{
	//			for (int cols = 0; cols < a.Cols(); cols++)
	//			{
	//				c[rows, cols] *= b;
	//			}
	//		}

	//		return c;
	//	}

	//	public static Matrix operator * (double a, Matrix b)
	//	{
	//		return b * a;
	//	}

	//	public static Matrix operator * (Matrix a, Matrix b)
	//	{
	//		if(a.Cols() != b.Rows())
	//		{
	//			return null;
	//		}
	//		else
	//		{
	//			var c = new Matrix(a.Rows(), b.Cols());

	//			for (int row = 0; row < c.Rows(); row++)
	//			{
	//				for (int col = 0; col < c.Cols(); col++)
	//				{
	//					double sum = 0;

	//					for (int step = 0; step < a.Cols(); step++)
	//					{
	//						sum += a[row, step] * b[step, col];
	//					}

	//					c[row, col] = sum;
	//				}
	//			}
	//			return c;
	//		}
	//	}

	//	public Matrix SubM(int row, int col)
	//	{
	//		var c = new Matrix(Rows() - 1, Cols() - 1);

	//		for (int rows = 0; rows < c.Rows(); rows++)
	//		{
	//			for (int cols = 0; cols < c.Cols(); cols++)
	//			{
	//				var nrow = (rows >= row) ? rows + 1 : rows;
	//				var ncol = (cols >= col) ? cols + 1 : cols;

	//				c[rows, cols] = m[nrow, ncol];
	//			}
	//		}

	//		return c;
	//	}

	//	public Matrix Transpose()
	//	{
	//		var o = new Matrix(Cols(), Rows());

	//		for(int row = 0; row < o.Rows(); row++)
	//		{
	//			for (int col = 0; col < o.Cols(); col++)
	//			{
	//				o[row, col] = this[col, row];
	//			}
	//		}

	//		return o;
	//	}

	//	public double Det()
	//	{
	//		if (Rows() == Cols()) {
	//			if (Rows() == 2 && Cols() == 2)
	//			{
	//				return m[0, 0] * m[1, 1] - m[1, 0] * m[0, 1];
	//			}
	//			else
	//			{
	//				double sum = 0;
	//				for (int col = 0; col < Cols(); col++)
	//				{
	//					sum += m[0, col] * SubM(0, col).Det();
	//				}

	//				return sum;
	//			}
	//		}
	//		else
	//		{
	//			return double.NaN;
	//		}
	//	}		
	//}

}
