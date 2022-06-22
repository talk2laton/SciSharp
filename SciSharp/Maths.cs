using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace SciSharp
{
    public class Maths
    {
        public static Format format = Format.Short;

        #region Transcedentals

        /// <summary>
        /// Returns the Absolute value of a number
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Abs(Nmbr v) => v.Abs();

        /// <summary>
        /// Returns the Argument value of a number
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Arg(Nmbr v) => v.Arg();

        /// <summary>
        /// Returns the Sign value of a number
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Sign(Nmbr v) => v.Sign();

        /// <summary>
        /// Returns the larger of two specified numbers.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Nmbr Max(Nmbr a, Nmbr b) => a > b ? a : b;

        /// <summary>
        /// Returns the smaller of two numbers.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Nmbr Min(Nmbr a, Nmbr b) => a < b ? a : b;

        /// <summary>
        /// Returns the square root of a specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Sqrt(Nmbr v) => v.Sqrt();

        /// <summary>
        /// Returns the square of a specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Sqr(Nmbr v) => v.Sqr();

        /// <summary>
        /// Returns the sine value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Sin(Nmbr v) => v.Sin();

        /// <summary>
        /// Returns the cosine value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Cos(Nmbr v) => v.Cos();

        /// <summary>
        /// Returns the tangent value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Tan(Nmbr v) => v.Tan();

        /// <summary>
        /// Returns the hyperbolic sine value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Sinh(Nmbr v) => v.Sinh();

        /// <summary>
        /// Returns the hyperbolic cosine value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Cosh(Nmbr v) => v.Cosh();

        /// <summary>
        /// Returns the hyperbolic tangent value of the specified angle
        /// </summary>
        /// <param name="v">v: an angle, measured in radians</param>
        /// <returns></returns>
        public static Nmbr Tanh(Nmbr v) => v.Tanh();

        /// <summary>
        /// Returns the angle whose sine is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Asin(Nmbr v) => v.Asin();

        /// <summary>
        /// Returns the angle whose cosine is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Acos(Nmbr v) => v.Acos();

        /// <summary>
        /// Returns the angle whose tangent is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Atan(Nmbr v) => v.Atan();

        /// <summary>
        /// Returns the angle whose hyperbolic sine is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Asinh(Nmbr v) => v.Asinh();

        /// <summary>
        /// Returns the angle whose hyperbolic cosine is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Acosh(Nmbr v) => v.Acosh();

        /// <summary>
        /// Returns the angle whose hyperbolic tangent is the specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Atanh(Nmbr v) => v.Atanh();

        /// <summary>
        /// Returns e raised to the specified power.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Exp(Nmbr v) => v.Exp();

        /// <summary>
        /// Returns the logarithm of a specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Log(Nmbr v) => v.Log();

        /// <summary>
        /// Returns the base 2 logarithm of a specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Log2(Nmbr v) => v.Log()/Log(2.0);

        /// <summary>
        /// Returns the base 10 logarithm of a specified number.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Nmbr Log10(Nmbr v) => v.Log() / Log(10.0);

        /// <summary>
        /// Returns the logarithm of value (v) to the base of the specified number (b).
        /// </summary>
        /// <param name="v"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Nmbr Log(Nmbr v, Nmbr b) => v.Log() / b.Log();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Nmbr Pow(Nmbr a, Nmbr b) => a ^ b;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr LnGamma(Nmbr x) => SpecialFunctions.LnGamma(x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr Gamma(Nmbr x) => SpecialFunctions.Gamma(x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr Erf(Nmbr x) => SpecialFunctions.Erf(x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr Erfc(Nmbr x) => SpecialFunctions.Erfc(x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr Zeta(Nmbr x) => SpecialFunctions.Zeta(x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr BesselJ(int n, Nmbr x) => SpecialFunctions.BesselClass.Jn(n, x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr BesselY(int n, Nmbr x) => SpecialFunctions.BesselClass.Yn(n, x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr BesselI(int n, Nmbr x) => SpecialFunctions.BesselClass.In(n, x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Nmbr BesselK(int n, Nmbr x) => SpecialFunctions.BesselClass.Kn(n, x);


        #endregion


        public enum Format
        {
            Long,
            Short
        }
    }
}
