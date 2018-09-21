using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Media;
using System.Collections;
using System.Windows.Media;
using System.Windows.Shapes;

namespace Utils
{
	// useful functions and constants and such.

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
			for (int i = 0; i < x.Dims; i++)
			{
				sum += Math.Exp(x[i]);
			}
			for (int i = 0; i < x.Dims; i++)
			{
				s.Add(Math.Exp(x[i] / sum));
			}
			return s;
		}

		public static Vec StiffMax(Vec x)
		{
			Vec s = new Vec();
			double sum = 0;
			for (int i = 0; i < x.Dims; i++)
			{
				sum += x[i];
			}
			for (int i = 0; i < x.Dims; i++)
			{
				s.Add(x[i] / sum);
			}
			return s;
		}

		public static int Choose(Vec x)
		{
			double c = rand.NextDouble();
			for (int i = 0; i < x.Dims; i++)
			{
				if (c < x[i])
				{
					return i;
				}
				c -= x[i];
			}
			return x.Dims;
		}

		public static int HotPick(Vec x)
		{
			double c = rand.NextDouble();
			int d = 0;
			while(c - x[d] > 0)
			{
				d++;
			}
			return d;
		}

		public static Vec Center(Vec a, Vec b, Vec c)
		{
			Matrix m = new Matrix();
			m.Add(new Vec() {   1,     1,    1,  1 });
			m.Add(new Vec() { a * a, a[0], a[1], 1 });
			m.Add(new Vec() { b * b, b[0], b[1], 1 });
			m.Add(new Vec() { c * c, c[0], c[1], 1 });
			Vec d = new Vec();
			d.Add((m.MinorMatrix(0, 1).Det / m.MinorMatrix(0, 0).Det) / 2);
			d.Add((m.MinorMatrix(0, 2).Det / m.MinorMatrix(0, 0).Det) / -2);
			return d;
		}

		public static Vec Center(Vec[] a)
		{
			return Center(a[0], a[1], a[2]);
		}

		public static bool IsInCirc(Vec[] a)
		{
			Matrix m = new Matrix();
			for (int i = 0; i < 4; i++)
			{
				m.Add(new Vec() { a[i][0], a[i][1], a[i] * a[i], 1 });
			}
			return m.Det < 0;
		}

		public static bool IsInPoly(System.Drawing.PointF[] polygon, System.Drawing.PointF testPoint)
		{
			bool result = false;
			int j = polygon.Count() - 1;
			for (int i = 0; i < polygon.Count(); i++)
			{
				if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
				{
					if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
					{
						result = !result;
					}
				}
				j = i;
			}
			return result;
		}

		public static bool IsInPoly(Point[] polygon, Point testPoint)
		{
			bool result = false;
			int j = polygon.Count() - 1;
			for (int i = 0; i < polygon.Count(); i++)
			{
				if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
				{
					if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
					{
						result = !result;
					}
				}
				j = i;
			}
			return result;
		}
		
		public static bool IsInPoly(Polygon p, Point testPoint)
		{
			bool result = false;

			var polygon = p.Points;

			int j = polygon.Count() - 1;
			for (int i = 0; i < polygon.Count(); i++)
			{
				if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
				{
					if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
					{
						result = !result;
					}
				}
				j = i;
			}
			return result;
		}

		public static bool IsInPoly(Polygon p, Vec testVector)
		{
			bool result = false;

			Point testPoint = new Point(testVector[0], testVector[1]);

			var polygon = p.Points;

			int j = polygon.Count() - 1;
			for (int i = 0; i < polygon.Count(); i++)
			{
				if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
				{
					if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
					{
						result = !result;
					}
				}
				j = i;
			}
			return result;
		}
	}

	// DataMatrix stuffs
	[Serializable]
	public class DataMatrix : IEnumerable
	{
		private List<dynamic> Data = new List<dynamic>();
		private List<double> Probabilies = new List<double>();



		public IEnumerator GetEnumerator()
		{
			yield return Data;
		}

		public void Add(params dynamic[] a)
		{
			Data.Add(a[0]);
			Probabilies.Add(a[1]);
		}

		

		public void Normalize()
		{
			double sum = 0;
			foreach(double d in Probabilies)
			{
				sum += d;
			}
			for (int p = 0; p < Probabilies.Count(); p++)
			{
				Probabilies[p] = Probabilies[p] / sum;
			}
		}

		public void InverseNormalize()
		{
			double sum = 0;
			foreach(double d in Probabilies)
			{
				sum += d;
			}
			for (int p = 0; p < Probabilies.Count(); p++)
			{
				Probabilies[p] = (sum - Probabilies[p]) / sum;
			}
		}

		public void SoftMax()
		{
			for (int p = 0; p < Probabilies.Count(); p++)
			{
				Probabilies[p] = Math.Exp(Probabilies[p]);
			}
			Normalize();
		}

		public dynamic Pick()
		{
			double c = Consts.rand.NextDouble();
			int p = 0;
			while((c -= Probabilies[p]) > 0)
			{
				p++;
				if(p >= Probabilies.Count)
				{
					return null;
				}
			}
			return Data[p];
		}
	}

	// Tensor stuffs
	[Serializable]
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
	[Serializable]
	public class Vec : IEnumerable
	{
		// vector variables

		private List<double> _v = new List<double>();

		// Vector element retrieval

		public double x
		{
			get
			{
				try
				{
					return _v[0];
				}
				catch
				{
					return double.NaN;
				}
			}
			set
			{
				_v[0] = value;
			}
		}

		public double y
		{
			get
			{
				try
				{
					return _v[1];
				}
				catch
				{
					return double.NaN;
				}
			}
			set
			{
				_v[1] = value;
			}
		}

		public double z
		{
			get
			{
				try { return _v[2]; }
				catch { return double.NaN; }
			}
			set
			{
				_v[2] = value;
			}
		}

		public double w
		{
			get
			{
				try { return _v[3]; }
				catch { return double.NaN; }
			}
			set
			{
				_v[3] = value;
			}
		}

		public double v
		{
			get
			{
				try { return _v[4]; }
				catch { return double.NaN; }
			}
			set
			{
				_v[4] = value;
			}
		}

		public double u
		{
			get
			{
				try { return _v[5]; }
				catch { return double.NaN; }
			}
			set
			{
				_v[5] = value;
			}
		}

		public double t
		{
			get
			{
				try { return _v[6]; }
				catch { return double.NaN; }
			}
			set
			{
				_v[6] = value;
			}
		}
		
		// Vector enumeration and indexing

		public IEnumerator GetEnumerator()
		{
			yield return _v;
		}

		public void Add(params double[] a)
		{
			_v.Add(a[0]);
		}

		public void Add(double d)
		{
			_v.Add(d);
		}

		public double this[int dim]
		{
			get
			{
				return _v[dim];
			}
			set
			{
				_v[dim] = value;
			}
		}

		// Vector constructors

		public Vec(int dims)
		{
			for (int dim = 0; dim < dims; dim++)
			{
				_v.Add(0);
			}
		}

		public Vec()
		{

		}

		public Vec(double[] a)
		{
			_v = a.ToList();
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
			foreach (double d in _v)
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
				return _v.Count();
			}
		}

		public Vec Copy()
		{
			var c = new Vec();
			foreach (double d in _v)
			{
				c._v.Add(d);
			}
			return c;
		}

		public double Mag
		{
			get => Math.Sqrt(this * this);
		}

		public void RemoveDim(int dim)
		{
			_v.RemoveAt(dim);
		}

		public void RemoveDims(int start, int num)
		{
			_v.RemoveRange(start, num);
		}

		public void Concat(Vec a)
		{
			_v.AddRange(a._v);
		}

		// Vector operations

		public static Vec operator + (Vec a, Vec b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c._v.Add(a[dim] + b[dim]);
			}
			return c;
		}

		public static Vec operator - (Vec a, Vec b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c._v.Add(a[dim] - b[dim]);
			}
			return c;
		}

		public static Vec operator * (Vec a, double b)
		{
			var c = new Vec();
			for (int dim = 0; dim < a.Dims; dim++)
			{
				c._v.Add(a[dim] * b);
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
				c._v.Add(a[dim] / b);
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
			return new Complex(a._v.ToArray());
		}

		public static explicit operator Quaterneon (Vec a)
		{
			return new Quaterneon(a._v.ToArray());
		}

		public static explicit operator Polar(Vec a)
		{
			Polar p = new Polar();
			p.r = a.Mag;
			p.t = Math.Atan2(a[1], a[0]);
			return p;
		}

		public static explicit operator Polar3(Vec a)
		{
			Polar3 p = new Polar3();
			p.r = a.Mag;
			p.h = Math.Atan(a[1] / a[0]);
			p.a = Math.Acos(a[2] / a.Mag);
			return p;
		}
	}

	// Polar stuffs
	[Serializable]
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

		public static Polar operator * (Polar a, Polar b)
		{
			return new Polar(a.r * b.r, a.t * b.t);
		}

		public static Polar operator / (Polar a, Polar b)
		{
			return new Polar(a.r / b.r, a.t / b.t);
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
	[Serializable]
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
	[Serializable]
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
	[Serializable]
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
	[Serializable]
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
			c.k = a.r * b.k + a.i * b.j - a.j * b.i + a.k * b.r;
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
	
}
