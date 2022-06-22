using System;
using System.Collections.Generic;
using System.Text;

namespace SciSharp
{
    public struct Nmbr : IComparable<Nmbr>
    {
        internal double x, y;
        internal bool isreal;

        public double Real { get { return x; } }
        public double Imag { get { return y; } }

        /// <summary>
        /// Constructing a Nmbr from a real number
        /// </summary>
        /// <param name="r"></param>
        public Nmbr(double r)
        { x = r; y = 0; isreal = true; }

        /// <summary>
        /// Constructing a Nmbr from a real and an imaginary number
        /// </summary>
        /// <param name="r"></param>
        /// <param name="i"></param>
        public Nmbr(double r, double i)
        { x = r; y = i; isreal = i == 0; }

        /// <summary>
        /// Explicit operator for conversion from Nmbr to double
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static explicit operator double(Nmbr value) => value.x;

        /// <summary>
        /// Explicit operator for conversion from Nmbr to Cmplx
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static explicit operator Cmplx(Nmbr value) => new Cmplx(value.x, value.y);

        /// <summary>
        /// Implicit operator for conversion from double to Nmbr
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static implicit operator Nmbr(double value) => new Nmbr(value);

        /// <summary>
        /// Implicit operator for conversion from Cmplx to Nmbr
        /// </summary>
        /// <param name="value"></param>
        public static implicit operator Nmbr(Cmplx value) => new Nmbr(value.x, value.y);

        /// <summary>
        /// Overrides the default Object.ToString() Method
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            double ay, xy = Math.Max(Math.Abs(x), ay = Math.Abs(y));
            int Dp = Maths.format == Maths.Format.Short ? 4 : 15, F = 4+Dp, mantissa = 0;
            if (xy != 0) mantissa = Math.Abs((int)Math.Floor(Math.Log(xy)));
            string zeros = "0000000000000000000000", fmtx, fmty, ans, pad = zeros.Substring(0, Dp);
            if (isreal)
            {
                fmtx = "{0," + F + ":0." + pad + "}";
                if (mantissa >= 2) fmtx = "{0," + F + ":0." + pad + "e+00}";
                ans = string.Format(fmtx, x) ;
            }
            else
            {
                if (mantissa >= 2)
                {
                    fmtx = "{0," + F + ":0." + pad + "e+00}";
                    fmty = "{0," + (-F) + ":0." + pad + "e+00i}";
                }
                else
                {
                    fmtx = "{0," + F + ":0." + pad + "}";
                    fmty = "{0," + (-F) + ":0." + pad + "i}";
                }
                ans = string.Format(fmtx, x) + (y >= 0 ? " + " : " - ") + string.Format(fmty, ay);
            }
            return ans;
        }

        /// <summary>
        /// Printing Nmbr in an Matrix
        /// </summary>
        /// <returns></returns>
        public string ToString4Array()
        {
            double ay = Math.Abs(y); int Dp = Maths.format == Maths.Format.Short ? 4 : 15, F = 4 + Dp;
            string zeros = "0000000000000000000000", ans, pad = zeros.Substring(0, Dp), 
                fmtx = "{0," + F + ":0." + pad + "}", fmty = "{0," + (-F) + ":0." + pad + "i}";
            ans = string.Format(fmtx, x);
            if (!isreal) ans += (y >= 0 ? " + " : " - ") + string.Format(fmty, ay);
            return ans;
        }

        /// <summary>
        /// Overrides Object.Equals()
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            Nmbr Other = (Nmbr)obj;
            if (ReferenceEquals(Other, null))
                return false;
            Nmbr This = this;
            if (This.isreal && Other.isreal)
                return x.Equals(Other.x);
            else
                return ((Cmplx)This).Equals((Cmplx)Other);
        }

        /// <summary>
        /// Overrides GetHashCode()
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            Nmbr This = this;
            if (This.isreal)
                return x.GetHashCode();
            else
                return x.GetHashCode() + y.GetHashCode();
        }

        /// <summary>
        /// Overrides Object.CompareTo()
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Nmbr other)
        {
            if (isreal && other.isreal)
                return x.CompareTo(other.x);
            else
                return ((Cmplx)this).CompareTo((Cmplx)other);
        }

        #region Transcedental functions
        internal Nmbr Abs()
        {
            Nmbr ans;
            if(isreal) ans = Math.Abs(x);
            else ans = Cmplx.Abs(new Cmplx(x,y));
            return ans;
        }
        internal Nmbr Arg()
        {
            Nmbr ans = 0.0;
            if(!isreal) ans = Cmplx.Arg((Cmplx)this);
            return ans;
        }
        internal Nmbr Sign()
        {
            Nmbr ans;
            if(isreal) ans = Math.Sign(x);
            else ans = Cmplx.Unit((Cmplx)this);
            return ans;
        }
        internal Nmbr Unit() => Sign();
        internal Nmbr Conj() => ~this;
        internal Nmbr Round(int decpts)
        {
            Nmbr ans;
            if(isreal) ans = Math.Round(x, decpts);
            else ans = Cmplx.Round((Cmplx)this, decpts);
            return ans;
        }
        internal Nmbr Sin()
        {
            Nmbr ans;
            if (isreal) ans = Math.Sin(x);
            else ans = Cmplx.Sin((Cmplx)this);
            return ans;
        }
        internal Nmbr Cos()
        {
            Nmbr ans;
            if (isreal) ans = Math.Cos(x);
            else ans = Cmplx.Cos((Cmplx)this);
            return ans;
        }
        internal Nmbr Tan()
        {
            Nmbr ans;
            if (isreal) ans = Math.Tan(x);
            else ans = Cmplx.Tan((Cmplx)this);
            return ans;
        }
        internal Nmbr Sinh()
        {
            Nmbr ans;
            if (isreal) ans = Math.Sinh(x);
            else ans = Cmplx.Sinh((Cmplx)this);
            return ans;
        }
        internal Nmbr Cosh()
        {
            Nmbr ans;
            if (isreal) ans = Math.Cosh(x);
            else ans = Cmplx.Cosh((Cmplx)this);
            return ans;
        }
        internal Nmbr Tanh()
        {
            Nmbr ans;
            if (isreal) ans = Math.Tanh(x);
            else ans = Cmplx.Tanh((Cmplx)this);
            return ans;
        }
        internal Nmbr Sqr()
        {
            Nmbr ans;
            if (isreal) ans = x * x;
            else ans = this * this;
            return ans;
        }
        internal Nmbr Exp()
        {
            Nmbr ans;
            if (isreal)
                ans = Math.Exp(x);
            else
            {
                Cmplx c = Cmplx.Exp((Cmplx)this);
                ans = new Nmbr(c.x, c.y);
            }
            return ans;
        }
        internal Nmbr Asin()
        {
            Nmbr ans;
            if (isreal)
            {
                if (Math.Abs(x) <= 1) ans = Math.Sin(x);
                else ans = Cmplx.Asin(x);
            }
            else ans = Cmplx.Asin((Cmplx)this);
            return ans;
        }
        internal Nmbr Acos()
        {
            Nmbr ans;
            if (isreal)
            {
                if (Math.Abs(x) <= 1) ans = Math.Acos(x);
                else ans = Cmplx.Acos(x);
            }
            else ans = Cmplx.Acos((Cmplx)this);
            return ans;
        }
        internal Nmbr Atan()
        {
            Nmbr ans;
            if (isreal) ans = Math.Atan(x);
            else ans = Cmplx.Atan((Cmplx)this);
            return ans;
        }
        internal Nmbr Asinh()
        {
            Nmbr ans;
            if (isreal) ans = Math.Log(x + Math.Sqrt(x * x + 1));
            else ans = Cmplx.Asinh((Cmplx)this);
            return ans;
        }
        internal Nmbr Acosh()
        {
            Nmbr ans;
            if (isreal)
            {
                if (x == 1) ans = 0.0;
                else if (x > 1) ans = Math.Log(x + Math.Sqrt(x * x - 1));
                else ans = Cmplx.Acosh(x);
            }
            else ans = Cmplx.Acosh((Cmplx)this);
            return ans;
        }
        internal Nmbr Atanh()
        {
            Nmbr ans;
            if (isreal)
            {
                if (Math.Abs(x) < 1) ans = 0.5 * Math.Log((1 + x) / (1 - x));
                else ans = Cmplx.Atanh(x);
            }
            else ans = Cmplx.Atanh((Cmplx)this);
            return ans;
        }
        internal Nmbr Sqrt()
        {
            Nmbr ans;
            if (isreal)
            {
                if (x == 0.0) ans = 0.0;
                else if (x > 0.0) ans = Math.Sqrt(x);
                else ans = Cmplx.Sqrt(x);
            }
            else ans = Cmplx.Sqrt((Cmplx)this);
            return ans;
        }
        internal Nmbr Log()
        {
            Nmbr ans;
            if (isreal)
            {
                if (x == 0.0) ans = double.NegativeInfinity;
                else if (x > 0.0) ans = Math.Log(x);
                else ans = Cmplx.Log(x);
            }
            else ans = Cmplx.Log((Cmplx)this);
            return ans;
        }
        #endregion

        #region operators
        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Nmbr operator +(Nmbr c1, Nmbr c2)
        {
            Nmbr ans;
            if (c1.isreal && c2.isreal)
            { ans = c1.x + c2.x;  }
            else
            { ans = (Cmplx)c1 + (Cmplx)c2; }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Nmbr operator -(Nmbr c)
        {
            Nmbr ans;
            if (c.isreal)
            { ans = new Nmbr(-c.x);  }
            else
            { ans = new Nmbr(-c.x, -c.y); }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Nmbr operator ~(Nmbr c)
        {
            Nmbr ans;
            if (c.isreal)
            { ans = new Nmbr(c.x); }
            else
            { ans = new Nmbr(c.x, -c.y); }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Nmbr operator -(Nmbr c1, Nmbr c2) => c1 + (-c2);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Nmbr operator *(Nmbr c1, Nmbr c2)
        {
            Nmbr ans;
            if (c1.isreal && c2.isreal)
            { double c = c1.x * c2.x; ans = new Nmbr(c); }
            else
            { Cmplx c = (Cmplx)c1 * (Cmplx)c2; ans = new Nmbr(c.x, c.y); }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Nmbr operator /(Nmbr c1, Nmbr c2) 
        {
            Nmbr ans;
            if (c1.isreal && c2.isreal)
            { double c = c1.x / c2.x; ans = new Nmbr(c); }
            else
            { Cmplx c = (Cmplx)c1 / (Cmplx)c2; ans = new Nmbr(c.x, c.y); }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Nmbr operator ^(Nmbr c1, Nmbr c2)
        {
            Nmbr ans;
            if (c1.isreal && c2.isreal)
            {
                if (c1.x > 0 || (int)c2.x - c2.x == 0)
                { double c = Math.Pow(c1.x, c2.x); ans = new Nmbr(c); }
                else
                { Cmplx c = (Cmplx)c1 ^ (Cmplx)c2; ans = new Nmbr(c.x, c.y); }
            }
            else
            { Cmplx c = (Cmplx)c1 ^ (Cmplx)c2; ans = new Nmbr(c.x, c.y); }
            return ans;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator ==(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) == 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator !=(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) != 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator >(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) > 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator <(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) < 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator >=(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) >= 0;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator <=(Nmbr c1, Nmbr c2) => c1.CompareTo(c2) <= 0;
        #endregion
    }
}
