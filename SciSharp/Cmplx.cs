using System;

namespace SciSharp
{
    public struct Cmplx : IComparable<Cmplx>
    {
        #region private
        internal double x, y;
        double GetMagnitude()
        {
            double a = Math.Abs(x), b = Math.Abs(y), result = 0.0;
            if (a == 0.0) result = b;
            else if (b == 0.0) result = a;
            else if (a != 0 && b != 0)
            {
                result = a > b ? a * Math.Sqrt(1 + Math.Pow(b / a, 2)) :
                                 b * Math.Sqrt(1 + Math.Pow(a / b, 2));
            }
            return result;
        }
        double GetArgument()
        {
            return Math.Atan2(y, x);
        }
        static Cmplx Conj(Cmplx c)
        {
            return new Cmplx(c.x, -c.y);
        }
        static Cmplx Times(Cmplx c1, Cmplx c2)
        {
            double a = c1.x, b = c1.y, c = c2.x, d = c2.y;
            double ac = a * c, bd = b * d, apb = a + b, cpd = c + d;
            Cmplx result = new Cmplx(ac - bd, apb * cpd - ac - bd);
            return result;
        }
        static Cmplx Div(Cmplx c1, Cmplx c2)
        {
            double a = c1.x, b = c1.y, numr = 0, numi, c = c2.x, d = c2.y, denr;
            if (Math.Abs(c) >= Math.Abs(d))
            {
                numr = a + b * (d / c);
                numi = b - a * (d / c);
                denr = c + d * (d / c);
            }
            else
            {
                numr = a * (c / d) + b;
                numi = b * (c / d) - a;
                denr = c * (c / d) + d;
            }
            return new Cmplx(numr / denr, numi / denr);
        }
        static Cmplx Pow(Cmplx c1, Cmplx c2)
        {
            Cmplx ans = Exp(Times(Log(c1), c2));
            return ans;
        }
        #endregion

        #region internal
        internal double Real { get { return x; } }
        internal double Imaginary { get { return y; } }
        internal double Magnitude { get { return GetMagnitude(); } }
        internal double Argument { get { return GetArgument(); } }
        internal static Cmplx Max(Cmplx A, Cmplx B)
        {
            int I = A.CompareTo(B);
            return I > 0 ? A : B;
        }
        internal static Cmplx Min(Cmplx A, Cmplx B)
        {
            int I = A.CompareTo(B);
            return I < 0 ? A : B;
        }
        internal static double Abs(Cmplx C)
        {
            return C.GetMagnitude();
        }
        internal static bool IsNaN(Cmplx C)
        {
            return double.IsNaN(C.Magnitude);
        }
        internal static bool IsReal(Cmplx C)
        {
            double c = Math.Abs(C.y) / C.Magnitude;
            return c < 1e-16;
        }
        internal static bool IsImaginary(Cmplx C)
        {
            double c = Math.Abs(C.x) / C.Magnitude;
            return c < 1e-16;
        }
        internal static Cmplx Round(Cmplx c, int n)
        {
            return new Cmplx(Math.Round(c.x, n), Math.Round(c.y, n));
        }
        internal static Cmplx Unit(Cmplx C)
        {
            double absC = Abs(C);
            return new Cmplx(C.x / absC, C.y / absC);
        }
        internal static double Arg(Cmplx c)
        {
            return Math.Atan2(c.y, c.x);
        }
        internal static Cmplx Cart(double abs, double angle)
        {
            double s = 0, c = 0;
            if (angle != Math.PI) s = Math.Sin(angle);
            if (angle != Math.PI / 2) c = Math.Cos(angle);
            return new Cmplx(abs * c, abs * s);
        }
        internal static double[] Pol(Cmplx c)
        {
            double[] ans = new double[2];
            ans[0] = Abs(c);
            ans[1] = Arg(c);
            return ans;
        }
        internal static Cmplx[] Root(Cmplx c, int N)
        {
            double p = 2 * Math.PI / N, angle = Arg(c) / N, v = Math.Pow(Abs(c), 1.0 / N);
            Cmplx[] ans = new Cmplx[N];
            for (int n = 0; n < N; n++)
            {
                ans[n] = Cart(v, angle + p * n);
            }
            return ans;
        }
        internal static Cmplx Sqrt(Cmplx c)
        {
            double x = c.x, y = c.y, w = 0, ax, ay, x2 = x * x, y2 = y * y; Cmplx ans;
            if ((ax = Math.Abs(x)) >= (ay = Math.Abs(y)))
                w = Math.Sqrt(ax) * Math.Sqrt((1 + Math.Sqrt(1 + (y2 / x2))) / 2);
            else if (ax < ay)
                w = Math.Sqrt(ay) * Math.Sqrt(((ax / ay) + Math.Sqrt(1 + (x2 / y2))) / 2);
            if (w == 0) ans = 0;
            else
            {
                if (x >= 0) ans = new Cmplx(w, y / (2 * w));
                else
                {
                    if (y >= 0)
                        ans = new Cmplx(ay / (2 * w), w);
                    else
                        ans = new Cmplx(ay / (2 * w), -w);
                }
            }
            return ans;
        }
        internal static Cmplx Exp(Cmplx c)
        {
            return Cart(Math.Exp(c.x), c.y);
        }
        internal static Cmplx Log(Cmplx c)
        {
            double[] p = Pol(c);
            Cmplx ans = new Cmplx(Math.Log(p[0]), p[1]);
            return ans;
        }
        internal static Cmplx Log(Cmplx c1, Cmplx c2)
        {
            return Log(c1) / Log(c2);
        }
        internal static Cmplx Sin(Cmplx c)
        {
            Cmplx I = new Cmplx(0, 1);
            Cmplx ans = (Exp(I * c) - Exp(-I * c)) / (2 * I);
            return ans;
        }
        internal static Cmplx Cos(Cmplx c)
        {
            Cmplx I = new Cmplx(0, 1);
            Cmplx ans = (Exp(I * c) + Exp(-I * c)) / (2);
            return ans;
        }
        internal static Cmplx Tan(Cmplx c)
        {
            Cmplx ans = Sin(c) / Cos(c);
            return ans;
        }
        internal static Cmplx Sinh(Cmplx c)
        {
            Cmplx ans = (Exp(c) - Exp(-c)) / (2);
            return ans;
        }
        internal static Cmplx Cosh(Cmplx c)
        {
            Cmplx ans = (Exp(c) + Exp(-c)) / (2);
            return ans;
        }
        internal static Cmplx Tanh(Cmplx c)
        {
            Cmplx ans = Sinh(c) / Cosh(c);
            return ans;
        }
        internal static Cmplx Asin(Cmplx c)
        {
            Cmplx I = new Cmplx(0, 1);
            Cmplx ans = (-I) * Log(I * c + Sqrt(1 - (c ^ 2)));
            return ans;
        }
        internal static Cmplx Acos(Cmplx c)
        {
            Cmplx I = new Cmplx(0, 1);
            Cmplx ans = (-I) * Log(c + I * Sqrt(1 - (c ^ 2)));
            return ans;
        }
        internal static Cmplx Atan(Cmplx c)
        {
            Cmplx I = new Cmplx(0, 1);
            Cmplx ans = (0.5 * I) * Log((1 - I * c) / (1 + I * c));
            return ans;
        }
        internal static Cmplx Asinh(Cmplx c)
        {
            Cmplx ans = Log(c + Sqrt(c * c + 1));
            return ans;
        }
        internal static Cmplx Acosh(Cmplx c)
        {
            Cmplx ans = Log(c + Sqrt((c * c) - 1));
            return ans;
        }
        internal static Cmplx Atanh(Cmplx c)
        {
            Cmplx ans = 0.5 * Log((1 + c) / (1 - c));
            return ans;
        }
        #endregion

        #region public
        

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Cmplx(double real, double imaginary)
        {
            x = real; y = imaginary;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Cmplx(double real)
        {
            x = real; y = 0;
        }

        /// <summary>
        /// Explicit operator for conversion from complex to double
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static explicit operator double(Cmplx value) => value.x;

        /// <summary>
        /// Implicit operator for conversion from double to complex
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static implicit operator Cmplx(double value) => new Cmplx(value, 0);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public override bool Equals(object other)
        {
            if (this == (Cmplx)other) return true;
            if (other == null) return false;
            return true;
        }

        public override int GetHashCode()
        {
            return x.GetHashCode() + y.GetHashCode();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Cmplx other)
        {
            // Alphabetic sort if Cmplx is equal. [A to Z]
            if (Magnitude == other.Magnitude)
                return Argument.CompareTo(other.Argument);
            
            // Default to Complex sort. [High to low]
            return Magnitude.CompareTo(other.Magnitude);
        }
        #endregion

        #region operators
        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Cmplx operator +(Cmplx c1, Cmplx c2) => new Cmplx(c1.x + c2.x, c1.y + c2.y);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Cmplx operator -(Cmplx c) => new Cmplx(-c.x, -c.y);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Cmplx operator ~(Cmplx c) => Conj(c);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Cmplx operator -(Cmplx c1, Cmplx c2) => c1 + (-c2);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Cmplx operator *(Cmplx c1, Cmplx c2) => Times(c1, c2);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Cmplx operator /(Cmplx c1, Cmplx c2) => Div(c1, c2);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Cmplx operator ^(Cmplx c1, Cmplx c2) => Pow(c1, c2);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator ==(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) == 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator !=(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) != 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator >(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) > 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator <(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) < 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator >=(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) >= 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator <=(Cmplx c1, Cmplx c2) => c1.CompareTo(c2) <= 0;
        #endregion
    }

    /// <summary>
    /// Exception from the Complex Class
    /// </summary>
    internal class CException : Exception
    {
        /// <summary>
        /// Exception from the Matrix Class
        /// </summary>
        /// <param name="Message">Message to be displayed</param>
        internal CException(string Message)
            : base(Message)
        { }
    }
}
