using System;
using System.Collections.Generic;
using System.IO.Pipes;
using System.Text;

namespace SciSharp
{
    /// <summary>
    /// Handles all computations invloving special functions
    /// </summary>
    internal class SpecialFunctions
    {
        static double[] a = new double[101];
        static double[] Fvals = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                               362880, 3628800, 39916800, 479001600,
                               6227020800, 87178291200, 1307674368000,
                               20922789888000, 355687428096000, 6402373705728000};
        static double[] BeN = { 1, -0.5, 1.0 / 6, 0, -1.0 / 30, 0, 1.0 / 42, 0, -1.0 / 30,
                                  0, 5.0 / 66, 0, -691.0 / 2730, 0, 7.0 / 6, 0, -3617.0 / 510 ,
                              0, 43867.0/798, 0, -174611.0/330, 0, 854513.0/138, 0, -236364091.0/2730,
                              0, -236364091.0/2730, 0, -23749461029.0/870, 0, 8615841276005.0/14322};
        static double EulerNo = 0.57721566490153286060651209008240243104215933593992;
        static Nmbr Fs = 0, Fc = 0, Fx = 0;

        /// <summary>
        /// Returns the value ln[Γ(xx)] for xx>0.
        /// </summary>
        /// <param name="xx">value of argument for which LnGamma is to be computed</param>
        /// <returns>ln[Γ(xx)]</returns>
        internal static Nmbr LnGamma(Nmbr xx)
        {
            Nmbr x, y, tmp, ser;
            Nmbr[] cof = {76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
            int j;
            y = xx;
            x = xx;
            tmp = x + 5.5;
            tmp -= (x + 0.5) * Maths.Log(tmp);
            ser = 1.000000000190015;
            for (j = 0; j <= 5; j++) { y += 1;  ser += cof[j] / y; }
            return -tmp + Maths.Log(2.5066282746310005 * ser / x);
        }

        /// <summary>
        /// Returns the value Γ(x).
        /// </summary>
        /// <param name="x">value of argument for which Gamma is to be computed</param>
        /// <returns>Γ(x)</returns>
        internal static Nmbr Gamma(Nmbr x)
        {
            if(x.isreal)
                return ((int)x - x == 0) ? Factorial((int)x.x - 1) : Maths.Exp(LnGamma(x));
            else
                return Maths.Exp(LnGamma(x));
        }

        /// <summary>
        /// Returns the value Γ_p(a,x).
        /// </summary>
        /// <param name="a">value of argument for which Upper Incomplete Gamma is to be computed</param>
        /// <param name="x">start of integration</param>
        /// <returns>Γ_p(a,x)</returns>
        internal static Nmbr GammaP(Nmbr a, Nmbr x)
        {
            //Returns the incomplete gamma functionP(a, x).
            if (x < 0.0 || a <= 0.0) throw new Exception("Invalid arguments in routine gammp");
            if (x < (a + 1.0))
            { //Use the series representation.
                Nmbr serln = gser(a, x);
                return serln.x;
            }
            else
            {
                //Use the continued fraction representation
                Nmbr cfln = gcf(a, x);
                return 1.0 - cfln.x; // and take its complement.
            }
        }

        /// <summary>
        /// Returns the value Γ_q(a,x).
        /// </summary>
        /// <param name="a">value of argument for which Lower Incomplete Gamma is to be computed</param>
        /// <param name="x">end of integration</param>
        /// <returns>Γ_p(a,x)</returns>
        internal static Nmbr GammaQ(Nmbr a, Nmbr x)
        {
            //Returns the incomplete gamma functionQ(a, x).
            if (x < 0.0 || a <= 0.0) throw new Exception("Invalid arguments in routine gammq");
            if (x < (a + 1.0))
            { //Use the series representation
                Nmbr serln = gser(a, x);
                return 1 - serln.x;
            }
            else
            { //Use the continued fraction representation.
                Nmbr cfln = gcf(a, x);
                return cfln.x; // and take its complement.
            }
        }

        /// <summary>
        /// Returns the value Γ_p^-1(p, a).
        /// </summary>
        /// <param name="a">value of argument for which Upper Incomplete Gamma is to be computed</param>
        /// <param name="x">start of integration</param>
        /// <returns>Γ_p^-1(p, a)</returns>
        internal static Nmbr InvGammaP(Nmbr p, Nmbr a)
        {
            //Returns x such that Γ_p(a,x) = p for an argument p between 0 and 1.
            int j;
            Nmbr x, err, t, u, pp, lna1 = 0, gln, afac = 0, a1 = a - 1;
            const double EPS = 1e-8; //Accuracy is the square of EPS.
            gln = LnGamma(a);
            if (a <= 0.0) throw new SpecialFunctionsException("a must be pos in InvGammaP");
            if (p >= 1.0) return Maths.Max(100.0, a + 100.0 * Maths.Sqrt(a));
            if (p <= 0.0) return 0.0;
            if (a > 1.0)
            {
                //Initial guess based on reference [1].
                lna1 = Maths.Log(a1);
                afac = Maths.Exp(a1 * (lna1 - 1.0) - gln);
                pp = (p < 0.5) ? p : 1.0 - p;
                t = Maths.Sqrt(-2.0 * Maths.Log(pp));
                x = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
                if (p < 0.5) x = -x;
                x = Maths.Max(1e-3, a * Maths.Pow(1.0 - 1.0 / (9.0 * a) - x / (3.0 * Maths.Sqrt(a)), 3));
            }
            else
            {
                //Initial guess based on equations (6.2.8) and (6.2.9). 
                t = 1.0 - a * (0.253 + a * 0.12);
                if (p < t) x = Maths.Pow(p / t, 1.0 / a);
                else x = 1.0 - Maths.Log(1.0 - (p - t) / (1.0 - t));
            }
            for (j = 0; j < 12; j++)
            {
                if (x <= 0.0) return 0.0; //x too small to compute accurately.
                err = GammaP(a, x) - p;
                if (a > 1.0) t = afac * Maths.Exp(-(x - a1) + a1 * (Maths.Log(x) - lna1));
                else t = Maths.Exp(-x + a1 * Maths.Log(x) - gln);
                u = err / t;
                x -= (t = u / (1.0 - 0.5 * Maths.Min(1.0, u * ((a - 1.0) / x - 1)))); //Halley’s method.
                if (x <= 0.0) x = 0.5 * (x + t);// Halve old value if x tries to go negative.
                if (Maths.Abs(t) < EPS * x) break;
            }
            return x;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="t"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        internal static Nmbr HeavisideTheta(Nmbr t, double a = 0.0)
        {
            return t >= a ? 1.0 : 0.0;
        }

        /// <summary>
        /// Error Function
        /// </summary>
        /// <param name="x">value of argument for which the error function is to be computed</param>
        /// <returns>Erf[x]</returns>
        internal static Nmbr Erf(Nmbr v)
        { //Returns the error function erf(x).
            if (v.isreal)
                return v < 0.0 ? -GammaP(0.5, v * v) : GammaP(0.5, v * v);
            else 
            {
                double x = v.x, y = v.y, x2 = x * x, xy = x * y, ky, fk, gk, pi = Math.PI; int kk = 0;
                Func<int, Nmbr> fgk = new Func<int, Nmbr>(k =>
                {
                    ky = k * y;
                    fk = 2 * x * (1 - Math.Cos(2 * xy) * Math.Cosh(ky)) + 
                                 k * Math.Sin(2 * xy) * Math.Sinh(ky);
                    gk = 2 * x * Math.Sin(2 * xy) * Math.Cosh(ky) +
                                 k * Math.Cos(2 * xy) * Math.Sinh(ky);
                    return new Nmbr(fk, gk);
                });
                Nmbr k2 = kk*kk, expmk2_4 = Maths.Exp(-k2/4), expmx2 = Maths.Exp(-x2),
                    e = Erf(x).x + expmx2 / pi * expmk2_4 / (k2 + 4 * x2) * fgk(kk);
                while(expmk2_4 > 1e-32)
                {
                    kk++; k2 = kk * kk; expmk2_4 = Maths.Exp(-k2 / 4); expmx2 = Maths.Exp(-x2);
                    e += 2* expmx2 / pi * expmk2_4 / (k2 + 4 * x2) * fgk(kk);
                }
                return e;
            }
        }

        /// <summary>
        /// Error Function
        /// </summary>
        /// <param name="x">value of argument for which the error function is to be computed</param>
        /// <returns>Erf[x]</returns>
        internal static Nmbr InvErf(Nmbr x)
        { //Returns a value whose error function is z.
            Nmbr a = 0.147, the_sign_of_x = Maths.Sign(x), z;
            if (x == 0)
            {
                return 0;
            }
            else
            {
                var ln_1minus_x_sqrd = Maths.Log(1 - x * x);
                var ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
                var ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
                var ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2 / (Math.PI * a));
                var first_sqrt = Maths.Sqrt((ln_etc_by2_plus2 * ln_etc_by2_plus2) - ln_1minusxx_by_a);
                var second_sqrt = Maths.Sqrt(first_sqrt - ln_etc_by2_plus2);
                return second_sqrt * the_sign_of_x;
            }
        }

        internal static Nmbr Laguerre(int n, Nmbr x)
        {
            Nmbr L0 = 1, L1 = -x + 1, L2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return L0;
            else if (n == 1) return L1;
            else
            {
                while (i < n)
                {
                    L2 = ((2.0 * i + 1.0 - x) * L1 - i * L0) / (i + 1);
                    L0 = L1; L1 = L2; i++;
                }
                return L2;
            }
        }

        internal static Nmbr Hermite(int n, Nmbr x)
        {
            Nmbr H0 = 1, H1 = 2 * x, H2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return H0;
            else if (n == 1) return H1;
            else
            {
                while (i < n)
                {
                    H2 = 2 * x * H1 - 2 * i * H0;
                    H0 = H1; H1 = H2; i++;
                }
                return H2;
            }
        }

        internal static Nmbr ChebyshevT(int n, Nmbr x)
        {
            Nmbr T0 = 1, T1 = x, T2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return T0;
            else if (n == 1) return T1;
            else
            {
                while (i < n)
                {
                    T2 = 2 * x * T1 - T0;
                    T0 = T1; T1 = T2; i++;
                }
                return T2;
            }
        }

        internal static Nmbr ChebyshevU(int n, Nmbr x)
        {
            Nmbr U0 = 1, U1 = 2 * x, U2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return U0;
            else if (n == 1) return U1;
            else
            {
                while (i < n)
                {
                    U2 = 2 * x * U1 - U0;
                    U0 = U1; U1 = U2; i++;
                }
                return U2;
            }
        }

        internal static Nmbr Legendre(int n, Nmbr x)
        {
            Nmbr P0 = 1, P1 = x, P2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return P0;
            else if (n == 1) return P1;
            else
            {
                while (i < n)
                {
                    P2 = 2 * i * P1 - P0 - (x * P1 - P0) / (i + 1);
                    P0 = P1; P1 = P2; i++;
                }
                return P2;
            }
        }

        /// <summary>
        /// Complemetary Error Function
        /// </summary>
        /// <param name="x">value of argument for which the complemetary error function is to be computed</param>
        /// <returns>Erfc[x]</returns>
        internal static Nmbr Erfc(Nmbr v)
        { //Returns the complementary error function erfc(x).
            if(v.isreal)
                return v < 0.0 ? 1.0 + GammaP(0.5, v * v) : GammaQ(0.5, v * v);
            else return 1 - Erf(v);
        }

        /// <summary>
        /// Zeta Function
        /// </summary>
        /// <param name="x">value of argument for which the zeta function is to be computed</param>
        /// <returns>Z[x]</returns>
        internal static Nmbr Zeta(Nmbr x)
        {
            Nmbr s = 1, sum = 1; int i = 1;
            if (x < 1)
            {
                return 1 - Zeta(2 - x);
            }
            else
            {
                while (s > 1e-20)
                {
                    i++;
                    s = Maths.Pow(i, -x);
                    sum += s;
                }
                return sum;
            }
        }

        /// <summary>
        /// x!
        /// </summary>
        /// <param name="x">value of argument for which factorial is to be computed</param>
        /// <returns>x!</returns>
        internal static Nmbr Factorial(int x)
        {
            return x < 19 ? Fvals[x] : Maths.Exp(LnGamma(x + 1));
        }

        /// <summary>
        /// Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B[a,b]</returns>
        internal static Nmbr Beta(Nmbr a, Nmbr b)
        {
            return Gamma(a) * Gamma(b) / Gamma(a + b);
        }

        /// <summary>
        /// Incomplete Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B[a, b, x]</returns>
        internal static Nmbr IncBeta(Nmbr a, Nmbr b, Nmbr x)
        {
            Nmbr bt; int SWITCH = 3000;
            if (a <= 0.0 || b <= 0.0) throw new SpecialFunctionsException("Bad a or b in routine betai");
            if (x < 0.0 || x > 1.0) throw new SpecialFunctionsException("Bad x in routine betai");
            if (x == 0.0 || x == 1.0) return x;
            if (a > SWITCH && b > SWITCH) return betaiapprox(a, b, x);
            bt = Maths.Exp(LnGamma(a + b) - LnGamma(a) - LnGamma(b) + a * Maths.Log(x) + b * Maths.Log(1.0 - x));
            if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
            else return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
        }

        /// <summary>
        /// Inverse Incomplete Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B^-1[p, a, b]</returns>
        internal static Nmbr InvIncBeta(Nmbr p, Nmbr a, Nmbr b)
        {
            return invbetai(p, a, b);
        }

        /// <summary>
        /// Exponential Integral(It is defined as a  definite integral of the ratio between an exponential function and powers of its argument.)
        /// </summary>
        /// <param name="n">power of the exponential argument</param>
        /// <param name="x">lower limit of integration</param>
        /// <returns>ExpInt_n[x]</returns>
        internal static Nmbr ExpInt(int n, Nmbr x)
        {
            int MAXIT = 100; //Maximum allowed number of iterations.
            double EULER = 0.5772156649;// Euler’s constantγ.
            double FPMIN = 1.0e-30; // Close to smallest representable floating-point number.
            double EPS = 1.0e-7; // Desired relative error, not smaller than the machine pre-cision.
            int i, ii, nm1;
            Nmbr a, b, c, d, del, fact, h, psi, ans;
            nm1 = n - 1;
            if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
                throw new Exception("bad arguments in expint");
            else
            {
                if (n == 0) ans = Maths.Exp(-x) / x; //Special case.
                else
                {
                    if (x == 0.0) ans = 1.0 / nm1; // Another special case.
                    else
                    {
                        if (x > 1.0)
                        { //Lentz’s algorithm (§5.2).
                            b = x + n;
                            c = 1.0 / FPMIN;
                            d = 1.0 / b;
                            h = d;
                            for (i = 1; i <= MAXIT; i++)
                            {
                                a = -i * (nm1 + i);
                                b += 2.0;
                                d = 1.0 / (a * d + b); //Denominators cannot be zero.
                                c = b + a / c;
                                del = c * d;
                                h *= del;
                                if (Maths.Abs(del - 1.0) < EPS)
                                {
                                    ans = h * Maths.Exp(-x);
                                    return ans;
                                }
                            }
                            throw new Exception("continued fraction failed in expint");
                        }
                        else
                        { //Evaluate series.
                            ans = (nm1 != 0 ? 1.0 / nm1 : -Maths.Log(x) - EULER); //Setfirstterm.
                            fact = 1.0;
                            for (i = 1; i <= MAXIT; i++)
                            {
                                fact *= -x / i;
                                if (i != nm1) del = -fact / (i - nm1);
                                else
                                {
                                    psi = -EULER; //Computeψ(n).
                                    for (ii = 1; ii <= nm1; ii++) psi += 1.0 / ii;
                                    del = fact * (-Maths.Log(x) + psi);
                                }
                                ans += del;
                                if (Maths.Abs(del) < Maths.Abs(ans) * EPS) return ans;
                            }
                            throw new Exception("series failed in expint");
                        }
                    }
                }
            }
            return ans;
        }

        /// <summary>
        /// Exponential Integral(It is defined as a  definite integral of the ratio between an exponential function and its argument.)
        /// </summary>
        /// <param name="x">Negative of lower limit of integration</param>
        /// <returns>EI[x] = -ExpInt[-x]</returns>
        internal static Nmbr EI(Nmbr x)
        {
            int MAXIT = 100; //Maximum allowed number of iterations.
            double EULER = 0.5772156649;// Euler’s constantγ.
            double FPMIN = 1.0e-30; // Close to smallest representable floating-point number.
            double EPS = 1.0e-8; // Desired relative error, not smaller than the machine pre-cision.
            int k;
            Nmbr fact, prev, sum, term;
            if (x <= 0.0) throw new Exception("Bad argument in ei");
            if (x < FPMIN) return Maths.Log(x) + EULER; //Special case: avoid failure of convergence test because of underflow. 
            if (x <= -Maths.Log(EPS))
            {
                sum = 0.0; // Use power series.
                fact = 1.0;
                for (k = 1; k <= MAXIT; k++)
                {
                    fact *= x / k;
                    term = fact / k;
                    sum += term;
                    if (term < EPS * sum) break;
                }
                if (k > MAXIT) throw new Exception("Series failed in ei");
                return sum + Maths.Log(x) + EULER;
            }
            else
            {// Use asymptotic series.
                sum = 0.0; //Start with second term.
                term = 1.0;
                for (k = 1; k <= MAXIT; k++)
                {
                    prev = term;
                    term *= k / x;
                    if (term < EPS) break;
                    //Since final sum is greater than one,termitself approximates the relative error.
                    if (term < prev) sum += term; //Still converging: add new term.
                    else
                    {
                        sum -= prev; //Diverging: subtract previous term and exit. 
                        break;
                    }
                }
                return Maths.Exp(x) * (1.0 + sum) / x;
            }
        }

        /// <summary>
        /// Sine Integral(It is defined as a definite integral of the ratio between an sine function and its argument.)
        /// </summary>
        /// <param name="x">Lower limit of integration</param>
        /// <returns>SI[x]</returns>
        internal static Nmbr SI(Nmbr x)
        {
            return cisi(x).x;
        }

        /// <summary>
        /// Cosine Integral(It is defined as a definite integral of the ratio between an cosine function and its argument.)
        /// </summary>
        /// <param name="x">Lower limit of integration</param>
        /// <returns>CI[x]</returns>
        internal static Nmbr CI(Nmbr x)
        {
            return cisi(x).y;
        }

        /// <summary>
        /// Fresnel Integral (It is defined as a definite integral of sine of square of its argument.)
        /// </summary>
        /// <param name="x">Upper limit of the integration</param>
        /// <returns>S[x]</returns>
        internal static Nmbr FresnelS(Nmbr x)
        {
            return frenel(x).x;
        }

        /// <summary>
        /// Fresnel Integral (It is defined as a definite integral of cosine of square of its argument.)
        /// </summary>
        /// <param name="x">Upper limit of the integration</param>
        /// <returns>C[x]</returns>
        internal static Nmbr FresnelC(Nmbr x)
        {
            return frenel(x).y;
        }

        internal static Nmbr EllipticF(Nmbr phi, Nmbr k)
        {
            Nmbr s = Maths.Sin(phi);
            return s * rf(Maths.Sqr(Maths.Cos(phi)), (1.0 - s * k) * (1.0 + s * k), 1.0);
        }

        internal static Nmbr EllipticK(Nmbr k)
        {
            Nmbr phi = Math.PI / 2;
            Nmbr s = Maths.Sin(phi);
            return s * rf(Maths.Sqr(Maths.Cos(phi)), (1.0 - s * k) * (1.0 + s * k), 1.0);
        }

        internal static Nmbr EllipticE(Nmbr phi, Nmbr k)
        {
            Nmbr s = Maths.Sin(phi), cc = Maths.Sqr(Maths.Cos(phi)), q = (1 - s * k) * (1 + s * k);
            return s * (rf(cc, q, 1.0) - Maths.Sqr(s * k) * rd(cc, q, 1.0) / 3.0);
        }

        internal static Nmbr EllipticE(Nmbr k)
        {
            Nmbr phi = Math.PI / 2;
            Nmbr s = Maths.Sin(phi), cc = Maths.Sqr(Maths.Cos(phi)), q = (1 - s * k) * (1 + s * k);
            return s * (rf(cc, q, 1.0) - Maths.Sqr(s * k) * rd(cc, q, 1.0) / 3.0);
        }

        internal static Nmbr EllipticPi(Nmbr phi, Nmbr n, Nmbr k)
        {
            Nmbr cc, enss, q, s;
            s = Maths.Sin(phi);
            enss = n * s * s;
            cc = Maths.Sqr(Maths.Cos(phi));
            q = (1.0 - s * k) * (1.0 + s * k);
            return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
        }

        internal static Nmbr EllipticPi(Nmbr n, Nmbr k)
        {
            Nmbr cc, enss, q, s, phi = Math.PI / 2;
            s = Maths.Sin(phi);
            enss = n * s * s;
            cc = Maths.Sqr(Maths.Cos(phi));
            q = (1.0 - s * k) * (1.0 + s * k);
            return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
        }

        #region private

        static Nmbr rf(Nmbr x, Nmbr y, Nmbr z)
        {
            //Computes Carlson’s elliptic integral of the first kind, RF(x,y,z)
            // x, y,and z must be non-negative, and at most one can be zero.
            Nmbr ERRTOL = 0.0025, THIRD = 1.0 / 3.0, C1 = 1.0 / 24.0, C2 = 0.1, C3 = 3.0 / 44.0, C4 = 1.0 / 14.0,
                TINY = 5.0 * 2.22507e-308, BIG = 0.2 * 1.79769e+308, alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty,
                sqrtz, xt, yt, zt;
            if (Maths.Min(Maths.Min(x, y), z) < 0.0 || Maths.Min(Maths.Min(x + y, x + z), y + z) < TINY || Maths.Max(Maths.Max(x, y), z) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rf");
            xt = x; yt = y; zt = z;
            do
            {
                sqrtx = Maths.Sqrt(xt);
                sqrty = Maths.Sqrt(yt);
                sqrtz = Maths.Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                ave = THIRD * (xt + yt + zt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
            } while (Maths.Max(Maths.Max(Maths.Abs(delx), Maths.Abs(dely)), Maths.Abs(delz)) > ERRTOL);
            e2 = delx * dely - delz * delz;
            e3 = delx * dely * delz;
            return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / Maths.Sqrt(ave);
        }

        static Nmbr rd(Nmbr x, Nmbr y, Nmbr z)
        {
            //Computes Carlson’s elliptic integral of the second kind, RD(x,y,z)
            // x and y must be non-negative, and at most one can be zero. z must be positive.
            Nmbr ERRTOL = 0.0015, C1 = 3.0 / 14.0, C2 = 1.0 / 6.0, C3 = 9.0 / 22.0, C4 = 3.0 / 26.0,
                C5 = 0.25 * C3, C6 = 1.5 * C4, TINY = 2.0 * Maths.Pow(1.79769e+308, -2.0 / 3), 
                BIG = 0.2 * ERRTOL * Maths.Pow(2.22507e-308, -2.0 / 3), alamb, ave, delx, dely, delz, 
                ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;
            if (Maths.Min(x, y) < 0.0 || Maths.Min(x + y, z) < TINY || Maths.Max(Maths.Max(x, y), z) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rd");
            xt = x; yt = y; zt = z; sum = 0.0; fac = 1.0;
            do
            {
                sqrtx = Maths.Sqrt(xt);
                sqrty = Maths.Sqrt(yt);
                sqrtz = Maths.Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                sum += fac / (sqrtz * (zt + alamb));
                fac = 0.25 * fac;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                ave = 0.2 * (xt + yt + 3.0 * zt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
            } while (Maths.Max(Maths.Max(Maths.Abs(delx), Maths.Abs(dely)), Maths.Abs(delz)) > ERRTOL);
            ea = delx * dely;
            eb = delz * delz;
            ec = ea - eb;
            ed = ea - 6.0 * eb;
            ee = ed + ec + ec;
            return 3.0 * sum + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee)
            + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea))) / (ave * Maths.Sqrt(ave));
        }

        static Nmbr rc(Nmbr x, Nmbr y)
        {
            //Computes Carlson’s degenerate elliptic integral, RC(x,y)
            // x must be non-negative, and y must be zero.
            // if y < 0 the cauchy principal is returned
            Nmbr ERRTOL = 0.0012, THIRD = 1.0 / 3.0, C1 = 0.3, C2 = 1.0 / 7.0, C3 = 0.375, C4 = 9.0 / 22.0,
                TINY = 5.0 * 2.22507e-308, BIG = 0.2 * 1.79769e+308,
                COMP1 = 2.236 / Maths.Sqrt(TINY), COMP2 = Maths.Sqr(TINY * BIG) / 25.0,
                     alamb, ave, s, w, xt, yt;
            if (x < 0.0 || y == 0.0 || (x + Maths.Abs(y)) < TINY || (x + Maths.Abs(y)) > BIG ||
                (y < -COMP1 && x > 0.0 && x < COMP2))
                throw new SpecialFunctionsException("invalid arguments in rd");
            xt = x; yt = y;
            if (y > 0.0)
            {
                xt = x; yt = y; w = 1.0;
            }
            else
            {
                xt = x - y; yt = -y; w = Maths.Sqrt(x) / Maths.Sqrt(xt);
            }
            do
            {
                alamb = 2.0 * Maths.Sqrt(xt) * Maths.Sqrt(yt) + yt;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                ave = THIRD * (xt + yt + yt);
                s = (yt - ave) / ave;
            } while (Maths.Abs(s) > ERRTOL);
            return w * (1.0 + s * s * (C1 + s * (C2 + s * (C3 + s * C4)))) / Maths.Sqrt(ave);
        }

        static Nmbr rj(Nmbr x, Nmbr y, Nmbr z, Nmbr p)
        {
            //Computes Carlson’s elliptic integral of the second kind, RJ(x,y,z,p)
            // x, y,and z must be nonnegative, and at most one can be zero.p must be nonzero. 
            //Ifp<0, the Cauchy principal value is returned.
            Nmbr ERRTOL = 0.0015, C1 = 3.0 / 14.0, C2 = 1.0 / 6.0, C3 = 9.0 / 22.0, C4 = 3.0 / 26.0,
                          C5 = 0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = C3 + C3, 
                          TINY = Maths.Pow(5 * 2.22507e-308, 1.0 / 3), BIG = 0.3 * Maths.Pow(1.79769e+308, 1.0 / 3),
                         a = 0, alamb, alpha, ans, ave, b = 0, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, 
                         pt, rcx = 0, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt, zt;
            if (Maths.Min(Maths.Min(x, y), z) < 0.0 || Maths.Min(Maths.Min(x + y, x + z), Maths.Min(y + z, Maths.Abs(p))) < TINY
                || Maths.Max(Maths.Max(x, y), Maths.Max(z, Maths.Abs(p))) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rj");
            sum = 0.0; fac = 1.0;
            if (p > 0)
            {
                xt = x; yt = y; zt = z; pt = p;
            }
            else
            {
                xt = Maths.Min(Maths.Min(x, y), z);
                zt = Maths.Max(Maths.Max(x, y), z);
                yt = x + y + z - xt - zt;
                a = 1.0 / (yt - p);
                b = a * (zt - yt) * (yt - xt);
                pt = yt + b;
                rho = xt * zt / yt;
                tau = p * pt / yt;
                rcx = rc(rho, tau);
            }
            do
            {
                sqrtx = Maths.Sqrt(xt);
                sqrty = Maths.Sqrt(yt);
                sqrtz = Maths.Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                alpha = Maths.Sqr(pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
                beta = pt * Maths.Sqr(pt + alamb);
                sum += fac * rc(alpha, beta);
                fac = 0.25 * fac;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                pt = 0.25 * (pt + alamb);
                ave = 0.2 * (xt + yt + zt + 2 * pt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
                delp = (ave - pt) / ave;
            } while (Maths.Max(Maths.Max(Maths.Abs(delx), Maths.Abs(dely)), Maths.Max(Maths.Abs(delz), Maths.Abs(delp))) > ERRTOL);
            ea = delx * (dely + delz) + dely * delz;
            eb = delx * dely * delz;
            ec = delp * delp;
            ed = ea - 3.0 * ec;
            ee = eb + 2.0 * delp * (ea - ec);
            ans = 3.0 * sum + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * ee) + eb * (C7 + delp * (-C8 + delp * C4))
            + delp * ea * (C2 - delp * C3) - C2 * delp * ec) / (ave * Maths.Sqrt(ave));
            if (p <= 0.0) ans = a * (b * ans + 3.0 * (rcx - rf(xt, yt, zt)));
            return ans;
        }

        static Nmbr betacf(Nmbr a, Nmbr b, Nmbr x)
        {
            //Evaluates continued fraction for incomplete beta function by 
            //modified Lentz’s method (5.2). User should not call directly.
            int m, m2;
            Nmbr aa, c, d, del, h, qab, qam, qap, EPS = 1e-52, FPMIN = 2.22507e-308 / EPS;
            qab = a + b; //Theseq’s will be used in factors that occur in the coefficients (6.4.6). 
            qap = a + 1.0;
            qam = a - 1.0;
            c = 1.0; //First step of Lentz’s method.
            d = 1.0 - qab * x / qap;
            if (Maths.Abs(d) < FPMIN) d = FPMIN;
            d = 1.0 / d;
            h = d;
            for (m = 1; m < 10000; m++)
            {
                m2 = 2 * m;
                aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d; // One step (the even one) of the recur-rence. 
                if (Maths.Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c;
                if (Maths.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1.0 + aa * d; //Next step of the recurrence (the oddone). 
                if (Maths.Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c;
                if (Maths.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = d * c;
                h *= del;
                if (Maths.Abs(del - 1.0) <= EPS) break; //Arewedone?
            }
            return h;
        }

        static Nmbr betaiapprox(Nmbr a, Nmbr b, Nmbr x)
        {
            //Evaluates Incomplete beta by quadrature 
            // User should not call directly.
            int j, m;
            Nmbr xu, t, sum, ans, a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
            Nmbr lnmu = Maths.Log(mu), lnmuc = Maths.Log(1.0 - mu);
            t = Maths.Sqrt(a * b / (Maths.Sqr(a + b) * (a + b + 1.0)));
            if (x > a / (a + b))
            { //Set how far to integrate into the tail:
                if (x >= 1.0) return 1.0;
                xu = Maths.Min(1, Maths.Max(mu + 10 * t, x + 5.0 * t));
            }
            else
            {
                if (x <= 0.0) return 0.0;
                xu = Maths.Max(0, Maths.Min(mu - 10 * t, x - 5.0 * t));
            }
            sum = 0;
            double[] y, w;
            gauleg(18, out y, out w);
            for (j = 0; j < 18; j++)
            { //Gauss-Legendre.
                t = x + (xu - x) * y[j];
                sum += w[j] * Maths.Exp(a1 * (Maths.Log(t) - lnmu) + b1 * (Maths.Log(1 - t) - lnmuc));
            }
            ans = sum * (xu - x) * Maths.Exp(a1 * lnmu - LnGamma(a) + b1 * lnmuc - LnGamma(b) + LnGamma(a + b));
            return ans > 0 ? 1 - ans : -ans;
        }

        static Nmbr invbetai(Nmbr p, Nmbr a, Nmbr b)
        {
            //Inverse of incomplete beta function. Returns xsuch that I_x(a,b) = p 
            //for argument p between 0 and 1.
            const double EPS = 1e-8;
            Nmbr pp, t, u, err, x = 0, al, h, w, afac, a1 = a - 1, b1 = b - 1;
            int j;

            if (p <= 0) return 0;
            else if (p >= 1) return 1;
            else if (a >= 1 && b >= 1)
            {
                pp = (p < 0.5) ? p : 1 - p;
                t = Maths.Sqrt(-2 * Maths.Log(pp));
                x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
                if (p < 0.5) x = -x;
                al = (Maths.Sqr(x) - 3.0) / 6.0;
                h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
                w = (x * Maths.Sqrt(al + h) / h) - (1.0 / (2.0 * b - 1) - 1.0 / (2.0 * a - 1.0)) * (al + 5.0 / 6.0 - 2.0 / (3.0 * h));
                x = a / (a + b * Maths.Exp(2 * w));
            }
            else
            {
                Nmbr lna = Maths.Log(a / (a + b)), lnb = Maths.Log(b / (a + b));
                t = Maths.Exp(a * lna) / a;
                u = Maths.Exp(b * lnb) / b;
                w = t + u;
                if (p < t / w) x = Maths.Pow(a * w * p, 1.0 / a);
                else x = 1.0 - Maths.Pow(b * w * (1.0 - p), 1.0 / b);
            }
            afac = -LnGamma(a) - LnGamma(b) + LnGamma(a + b);
            for (j = 0; j < 10; j++)
            {
                if (x == 0.0 || x == 1.0) return x; //a or b too small for accurate calcu-lation. 
                err = IncBeta(a, b, x) - p;
                t = Maths.Exp(a1 * Maths.Log(x) + b1 * Maths.Log(1.0 - x) + afac);
                u = err / t; //Halley:
                x -= (t = u / (1.0 - 0.5 * Maths.Min(1.0, u * (a1 / x - b1 / (1.0 - x)))));
                if (x <= 0.0) x = 0.5 * (x + t); //Bisect if xtries to go neg or>1.
                if (x >= 1.0) x = 0.5 * (x + t + 1.0);
                if (Maths.Abs(t) < EPS * x && j > 0) break;
            }
            return x;
        }

        static void gauleg(int n, out double[] x, out double[] w)
        {
            // Given the lower and upper limits of integrationx1andx2,andgivenn, 
            // this routine returns arraysx[1..n]andw[1..n]of lengthn, 
            // containing the abscissas and weights of the Gauss-Legendren-point quadrature formula.
            double EPS = 1e-10;
            int m, j, i;
            double x1 = -1, x2 = 1, z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this rou-tine.
            x = new double[n]; w = new double[n];
            m = (n + 1) / 2; //The roots are symmetric in the interval, so
            //we only have to find half of them. 
            xm = 0.5 * (x1 + x2);
            xl = 0.5 * (x2 - x1);
            for (i = 0; i < m; i++)
            { //Loop over the desired roots.
                z = Math.Cos(Math.PI * (i + 1 - 0.25) / (n + 0.5));
                //Starting with the above approximation to theith root, we enter the main loop of
                //refinement by Newton’s method.
                do
                {
                    p1 = 1.0;
                    p2 = 0.0;
                    for (j = 1; j <= n; j++)
                    { //Loop up the recurrence relation to get the
                        //Legendre polynomial evaluated atz. 
                        p3 = p2;
                        p2 = p1;
                        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
                    }
                    //p1is now the desired Legendre polynomial. We next compute pp, its derivative,
                    //by a standard relation involving alsop2, the polynomial of one lower order.
                    pp = n * (z * p1 - p2) / (z * z - 1.0);
                    z1 = z;
                    z = z1 - p1 / pp; //Newton’s method.
                } while (Math.Abs(z - z1) > EPS);
                x[i] = xm - xl * z; // Scale the root to the desired interval,
                x[n - 1 - i] = xm + xl * z; // and put in its symmetric counterpart.
                w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp); // Compute the weight
                w[n - 1 - i] = w[i]; // and its symmetric counterpart.
            }
        }

        static Nmbr BesselJ0(Nmbr x)
        {
            Nmbr ax, xxx, xx, y, ans, ans1, ans2;
            ax = Maths.Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                       + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                        + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                     + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                       + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                        - y * 0.934945152e-7)));
                ans = Maths.Sqrt(0.636619772 / ax) * (Maths.Cos(xx) * ans1 - xxx * Maths.Sin(xx) * ans2);
            }
            return ans;
        }

        static Nmbr BesselJ1(Nmbr x)
        {
            Nmbr ax, xxx, xx, y, ans, ans1, ans2;
            ax = Maths.Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                       + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                       + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                      + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
                ans = Maths.Sqrt(0.636619772 / ax) * (Maths.Cos(xx) * ans1 - xxx * Maths.Sin(xx) * ans2);
                if (x.x < 0)ans = -ans;
            }
            return ans;
        }

        static Nmbr BesselY0(Nmbr x)
        {
            Nmbr xxx, xx, y, ans, ans1, ans2;
            if (x < 8.0)
            {
                y = x * x;
                ans1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
                       + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
                ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                       + y * (47447.26470 + y * (226.1030244 + y * 1.0))));
                ans = (ans1 / ans2) + 0.636619772 * BesselJ0(x) * Maths.Log(x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                        + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                        + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
                ans = Maths.Sqrt(0.636619772 / x) * (Maths.Sin(xx) * ans1 + xxx * Maths.Cos(xx) * ans2);
            }
            return ans;
        }

        static Nmbr BesselY1(Nmbr x)
        {
            Nmbr xxx, xx, y, ans, ans1, ans2;
            if (x < 8.0)
            {
                y = x * x;
                ans1 = x * (-0.4900604943e13 + y * (0.1275274390e13
                        + y * (-0.5153438139e11 + y * (0.7349264551e9
                        + y * (-0.4237922726e7 + y * 0.8511937935e4)))));
                ans2 = 0.2499580570e14 + y * (0.4244419664e12
                       + y * (0.3733650367e10 + y * (0.2245904002e8
                       + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
                ans = (ans1 / ans2) + 0.636619772 * (BesselJ1(x) * Maths.Log(x) - 1.0 / x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                       + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6
                       + y * 0.105787412e-6)));
                ans = Maths.Sqrt(0.636619772 / x) * (Maths.Sin(xx) * ans1 + xxx * Maths.Cos(xx) * ans2);
            }
            return ans;
        }

        static Nmbr BesselI0(Nmbr x)
        {
            Nmbr ax, ans, y; // Accumulate polynomials in double precision.
            if ((ax = Maths.Abs(x)) < 3.75) //Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                       + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
            }
            else
            {
                y = 3.75 / ax;
                ans = (Maths.Exp(ax) / Maths.Sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                     + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                     + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                     + y * 0.392377e-2))))))));
            }
            return ans;
        }

        static Nmbr BesselI1(Nmbr x)
        {
            Nmbr ax, ans, y;
            ax = Maths.Abs(x);
            if (ax < 3.75) // Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                    + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
            }
            else
            {
                y = 3.75 / ax;
                ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
                ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
                ans *= (Maths.Exp(ax) / Maths.Sqrt(ax));
            }
            if (x < 0) ans = -ans;
            return ans;
        }

        static Nmbr BesselK0(Nmbr x)
        {
            Nmbr y, ans; //Accumulate polynomials in double precision
            if (x <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (-Maths.Log(x / 2.0) * BesselI0(x)) + (-0.57721566 + y * (0.42278420
                    + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
                    + y * (0.10750e-3 + y * 0.74e-5))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Maths.Exp(-x) / Maths.Sqrt(x)) * (1.25331414 + y * (-0.7832358e-1
                    + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                    + y * (-0.251540e-2 + y * 0.53208e-3))))));
            }
            return ans;
        }

        static Nmbr BesselK1(Nmbr x)
        {
            Nmbr y, ans;
            if (x <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (Maths.Log(x / 2.0) * BesselI1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                    + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
                    + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Maths.Exp(-x) / Maths.Sqrt(x)) * (1.25331414 + y * (0.23498619
                    + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                    + y * (0.325614e-2 + y * (-0.68245e-3)))))));
            }
            return ans;
        }

        static Nmbr Chebev(Nmbr a, Nmbr b, Nmbr[] c, Nmbr x)
        {
            Nmbr d = 0.0, dd = 0.0, sv, y, y2;
            int j, m = c.Length;
            if ((x - a) * (x - b) > 0.0) throw new Exception("x not in range in routine chebev");
            y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a)); //Change of variable.
            for (j = m - 1; j >= 1; j--) // Clenshaw’s recurrence.
            {
                sv = d;
                d = y2 * d - dd + c[j];
                dd = sv;
            }
            return y * d - dd + 0.5 * c[0]; //Last step is different.
        }

        static void beschb(Nmbr x, out Nmbr gam1, out Nmbr gam2, out Nmbr gampl, out Nmbr gammi)
        {
            //Evaluates Γ1 and Γ2 by Chebyshev expansion for|x|≤1/2. Also returns 1/Γ(1 +x) and
            //1/Γ(1−x). If converting to double precision, set NUSE1=7, NUSE2=8.
            //double chebev(double a, double b, double[] c, float x);
            Nmbr xx;
            Nmbr[] c1 = {-1.142022680371168e0, 6.5165112670737e-3, 3.087090173086e-4,
                                  -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13};
            Nmbr[] c2 = {1.843740587300905e0,-7.68528408447867e-2,1.2719271366546e-3,
                              -4.9717367042e-6,-3.31261198e-8,2.423096e-10,-1.702e-13,-1.49e-15};
            xx = 8.0 * x * x - 1.0; // Multiply x by 2 to make range be −1 to 1,
            // and then apply transformation for eval-uating even Chebyshev series.
            gam1 = Chebev(-1.0, 1.0, c1, xx);
            gam2 = Chebev(-1.0, 1.0, c2, xx);
            gampl = gam2 - x * gam1;
            gammi = gam2 + x * (gam1);
        }

        static void besseljy(Nmbr x, Nmbr xnu, out Nmbr rj, out Nmbr ry, out Nmbr rjp, out Nmbr ryp)
        {
            Nmbr EPS = 1.0e-16, FPMIN = 1.0e-30, XMIN = 2.0, PI = Math.PI;
            int MAXIT = 10000;
            /*
             Returns the Bessel functionsrj=Jν, ry=Yν and their derivativesrjp=Jν, ryp=Yν,for
             positivexand forxnu=ν≥0. The relative accuracy is within one or two significant digits
             ofEPS, except near a zero of one of the functions, where EPScontrols its absolute accuracy.
             FPMINis a number close to the machine’s smallest floating-point number. All internal arithmetic
             is in double precision.
             */

            int i, isign, l, nl;
            Nmbr a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
            fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
            rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
            temp, w, x2, xi, xi2, xmu, xmu2;

            nl = x < XMIN ? (int)(xnu + 0.5) : Math.Max(0, (int)(xnu - x + 1.5));
            xmu = xnu - nl;
            xmu2 = xmu * xmu;
            xi = 1.0 / x;
            xi2 = 2.0 * xi;
            w = xi2 / PI;    // The Wronskian.
            isign = 1;      // Evaluate CF1 by modified Lentz’s method (§5.2).
            // isignkeeps track of sign changes in the de-nominator.
            h = Maths.Max(xnu * xi, FPMIN);
            b = xi2 * xnu;
            d = 0.0;
            c = h;
            for (i = 1; i <= MAXIT; i++)
            {
                b += xi2;
                d = b - d;
                if (Maths.Abs(d) < FPMIN) d = FPMIN;
                c = b - 1.0 / c;
                if (Maths.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = c * d;
                h = del * h;
                if (d < 0.0) isign = -isign;
                if (Maths.Abs(del - 1.0) < EPS) break;
            }

            if (i > MAXIT) throw new Exception("x too large in bessjy; try asymptotic expansion");
            rjl = isign * FPMIN; // Initialize Jν and Jν for downward recurrence.
            rjpl = h * rjl;
            rjl1 = rjl; //Store values for later rescaling.
            rjp1 = rjpl;
            fact = xnu * xi;
            for (l = nl; l >= 1; l--)
            {
                rjtemp = fact * rjl + rjpl;
                fact -= xi;
                rjpl = fact * rjtemp - rjl;
                rjl = rjtemp;
            }
            if (rjl == 0.0) rjl = EPS;
            f = rjpl / rjl; // Now have unnormalized Jµ and Jµ.
            if (x < XMIN) //Use series.
            {
                x2 = 0.5 * x;
                pimu = PI * xmu;
                fact = (Maths.Abs(pimu) < EPS ? 1.0 : pimu / Maths.Sin(pimu));
                d = -Maths.Log(x2);
                e = xmu * d;
                fact2 = (Maths.Abs(e) < EPS ? 1.0 : Maths.Sinh(e) / e);
                beschb(xmu, out gam1, out gam2, out gampl, out gammi); //Chebyshev evaluation ofΓ1andΓ2.
                ff = 2.0 / PI * fact * (gam1 * Maths.Cosh(e) + gam2 * fact2 * d); //f0.
                e = Maths.Exp(e);
                p = e / (gampl * PI); //p0.
                q = 1.0 / (e * PI * gammi); //q0.
                pimu2 = 0.5 * pimu;
                fact3 = (Maths.Abs(pimu2) < EPS ? 1.0 : Maths.Sin(pimu2) / pimu2);
                r = PI * pimu2 * fact3 * fact3;
                c = 1.0;
                d = -x2 * x2;
                sum = ff + r * q;
                sum1 = p;
                for (i = 1; i <= MAXIT; i++)
                {
                    ff = (i * ff + p + q) / (i * i - xmu2);
                    c *= (d / i);
                    p /= (i - xmu);
                    q /= (i + xmu);
                    del = c * (ff + r * q);
                    sum += del;
                    del1 = c * p - i * del;
                    sum1 += del1;
                    if (Maths.Abs(del) < (1.0 + Maths.Abs(sum)) * EPS) break;
                }
                if (i > MAXIT) throw new Exception("bessy series failed to converge");
                rymu = -sum;
                ry1 = -sum1 * xi2;
                rymup = xmu * xi * rymu - ry1;
                rjmu = w / (rymup - f * rymu); //Equation (6.7.13).
            }
            else // Evaluate CF2 by modified Lentz’s method (§5.2).
            {
                a = 0.25 - xmu2;
                p = -0.5 * xi;
                q = 1.0;
                br = 2.0 * x;
                bi = 2.0;
                fact = a * xi / (p * p + q * q);
                cr = br + q * fact;
                ci = bi + p * fact;
                den = br * br + bi * bi;
                dr = br / den;
                di = -bi / den;
                dlr = cr * dr - ci * di;
                dli = cr * di + ci * dr;
                temp = p * dlr - q * dli;
                q = p * dli + q * dlr;
                p = temp;
                for (i = 2; i <= MAXIT; i++)
                {
                    a += 2 * (i - 1);
                    bi += 2.0;
                    dr = a * dr + br;
                    di = a * di + bi;
                    if (Maths.Abs(dr) + Maths.Abs(di) < FPMIN) dr = FPMIN;
                    fact = a / (cr * cr + ci * ci);
                    cr = br + cr * fact;
                    ci = bi - ci * fact;
                    if (Maths.Abs(cr) + Maths.Abs(ci) < FPMIN) cr = FPMIN;
                    den = dr * dr + di * di;
                    dr /= den;
                    di /= -den;
                    dlr = cr * dr - ci * di;
                    dli = cr * di + ci * dr;
                    temp = p * dlr - q * dli;
                    q = p * dli + q * dlr;
                    p = temp;
                    if (Maths.Abs(dlr - 1.0) + Maths.Abs(dli) < EPS) break;
                }
                if (i > MAXIT) throw new Exception("cf2 failed in bessjy");
                gam = (p - f) / q; // Equations (6.7.6) – (6.7.10).
                rjmu = Maths.Sqrt(w / ((p - f) * gam + q));
                rjmu = Maths.Sign(rjmu) * Maths.Abs(rjl);
                rymu = rjmu * gam;
                rymup = rymu * (p + q / gam);
                ry1 = xmu * xi * rymu - rymup;
            }
            fact = rjmu / rjl;
            rj = rjl1 * fact; //Scale originalJν andJν
            rjp = rjp1 * fact;
            for (i = 1; i <= nl; i++) // Upward recurrence ofYν.
            {
                rytemp = (xmu + i) * xi2 * ry1 - rymu;
                rymu = ry1;
                ry1 = rytemp;
            }
            ry = rymu;
            ryp = xnu * xi * rymu - ry1;
        }

        static Nmbr gser(Nmbr a, Nmbr x)
        {
            const int ITMAX = 100;
            const double EPS = 3.0e-7;
            Nmbr sum, del, ap, gamser, gln;
            Nmbr ans; 
            int n;
            gln = LnGamma(a);
            ans = new Nmbr(0, gln.x);
            if (x <= 0.0)
            {
                if (x < 0.0) throw new Exception("x less than 0 in routine gser");
                gamser = 0.0;
                ans += gamser;
                return ans;
            }
            else
            {
                ap = a;
                del = sum = 1.0 / a;
                for (n = 1; n <= ITMAX; n++)
                {
                    ap += 1;
                    del *= x / ap;
                    sum += del;
                    if (Maths.Abs(del) < Maths.Abs(sum) * EPS)
                    {
                        gamser = sum * Maths.Exp(-x + a * Maths.Log(x) - gln);
                        ans += gamser;
                        return ans;
                    }
                }
                throw new Exception("a too large, ITMAX too small in routine gser");
            }
        }

        static Nmbr gcf(Nmbr a, Nmbr x)
        {
            const int ITMAX = 100;
            const double EPS = 3.0e-7, FPMIN = 1e-30;
            int i;
            Nmbr an, b, c, d, del, h, gln, gammcf, ans;
            gln = LnGamma(a); ans = new Nmbr(0, gln.x);
            b = x + 1.0 - a; //Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0=0.
            c = 1.0 / FPMIN;
            d = 1.0 / b;
            h = d;
            for (i = 1; i <= ITMAX; i++)
            { //Iterate to convergence.
                an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (Maths.Abs(d) < FPMIN) d = FPMIN;
                c = b + an / c;
                if (Maths.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = d * c;
                h *= del;
                if (Maths.Abs(del - 1.0) < EPS) break;
            }
            gammcf = Maths.Exp(-x + a * Maths.Log(x) - gln) * h;
            ans += gammcf;
            return ans;
        }

        static Nmbr cisi(Nmbr x)
        {
            // Computes the cosine and sine integrals Ci(x) and Si(x). The function Ci(x) is returned 
            // as the real part of the cs, and Si(x) as the imaginary part. Ci(0) is returned as a 
            // large negative number and no error message is generated. For x < 0 the routine returns 
            // Ci(−x) and you must supply the −iπ (which is the Si(x)) yourself.
            double EPS = 6.0e-8; // Relative error, or absolute error near a zero ofCi(x).
            double EULER = 0.57721566; //Euler’s constant γ.
            int MAXIT = 100; // Maximum number of iterations allowed.
            double PIBY2 = 1.5707963; // π/2.
            double FPMIN = 1.0e-30; //Close to smallest representable floating-point number.
            double TMIN = 2.0; //Dividing line between using the series and continued fraction. 
            int i, k;
            bool odd;
            Nmbr a, err, fact, sign, sum, sumc = 0, sums, t, term, si, ci;
            Nmbr h, b, c, d, del, cs; 
            t = Maths.Abs(x);
            if (t == 0.0) cs = -1.0 / FPMIN; // Special case.
            if (t > TMIN)
            { //Evaluate continued fraction by modified
                //Lentz’s method (§5.2). 
                b = new Cmplx(1.0, t.x);
                c = 1.0 / FPMIN;
                d = h = 1 / b;
                for (i = 2; i <= MAXIT; i++)
                {
                    a = -(i - 1) * (i - 1);
                    b += 2;
                    d = 1 / ((a * d) + b); //Denominators cannot be zero.
                    c = b + (a / c);
                    del = c * d;
                    h *= del;
                    if (Maths.Abs(del.x - 1.0) + Maths.Abs(del.y) < EPS) break;
                }
                if (i > MAXIT) throw new Exception("cf failed in cisi");
                h *= Cmplx.Cart(1, -t.x);
                cs = -~h + new Cmplx(0.0, PIBY2);
            }
            else
            { //Evaluate both series simultaneously.
                if (t < Math.Sqrt(FPMIN))
                { //Special case: avoid failure of convergence test because of underflow. sumc=0.0;
                    sums = t;
                }
                else
                {
                    sum = sums = sumc = 0.0;
                    sign = fact = 1.0;
                    odd = true;
                    for (k = 1; k <= MAXIT; k++)
                    {
                        fact *= t / k;
                        term = fact / k;
                        sum += sign * term;
                        err = term / Maths.Abs(sum);
                        if (odd)
                        {
                            sign = -sign;
                            sums = sum;
                            sum = sumc;
                        }
                        else
                        {
                            sumc = sum;
                            sum = sums;
                        }
                        if (err < EPS) break;
                        odd = !odd;
                    }
                    if (k > MAXIT) throw new Exception("maxits exceeded in cisi");
                }
                si = sums;
                ci = sumc + Maths.Log(t) + EULER;
                cs = new Cmplx(ci.x, si.x);
            }
            if (x < 0.0) { cs = ~cs; }
            return cs;
        }

        static Nmbr frenel(Nmbr x)
        {
            // Computes the Fresnel integralsS(x) andC(x) for all real x.
            Nmbr EPS = 6.0e-8; // Relative error, or absolute error near a zero ofCi(x).
            Nmbr EULER = 0.57721566; //Euler’s constant γ.
            int MAXIT = 100; // Maximum number of iterations allowed.
            Nmbr PIBY2 = 1.5707963; // π/2.
            Nmbr PI = 3.1415927; // π.
            Nmbr FPMIN = 1.0e-30; //Close to smallest representable floating-point number.
            Nmbr XMIN = 1.5; //Dividing line between using the series and continued fraction.
            int k, n;
            bool odd;
            Nmbr a, ax, fact, pix2, sign, sum, sumc, sums, term, test, s, c, b, cc, d, h, del, cs;
            ax = Maths.Abs(x);
            if (ax < Maths.Sqrt(FPMIN))
            {// Special case.
                cs = ax;
            }
            else if (ax < XMIN)
            {
                sum = sums = 0.0;
                sumc = ax;
                sign = 1.0;
                fact = PIBY2 * ax * ax;
                odd = true;
                term = ax;
                n = 3;
                for (k = 1; k <= MAXIT; k++)
                {
                    term *= fact / k;
                    sum += sign * term / n;
                    test = Maths.Abs(sum) * EPS;
                    if (odd)
                    {
                        sign = -sign;
                        sums = sum;
                        sum = sumc;
                    }
                    else
                    {
                        sumc = sum;
                        sum = sums;
                    }
                    if (term < test) break;
                    odd = !odd;
                    n += 2;
                }
                if (k > MAXIT) throw new Exception("series failed in frenel");
                cs = new Cmplx(sumc.x, sums.x);
            }
            else
            { //Evaluate continued fraction by modified Lentz’s method (§5.2). 
                pix2 = PI * ax * ax;
                b = new Cmplx(1.0, -pix2.x);
                cc = 1.0 / FPMIN;
                d = h = 1 / b;
                n = -1;
                for (k = 2; k <= MAXIT; k++)
                {
                    n += 2;
                    a = -n * (n + 1);
                    b += 4;
                    d = 1.0 / (a * d + b); //Denominators cannot be zero.
                    cc = b + a / cc;
                    del = cc * d;
                    h *= del;
                    if (Maths.Abs(del.x - 1.0) + Maths.Abs(del.y) < EPS) break;
                }
                if (k > MAXIT) throw new Exception("cf failed in frenel");
                h *= new Cmplx(ax.x, -ax.x);
                cs = new Cmplx(0.5, 0.5) * (1 - Cmplx.Cart(1, 0.5 * pix2.x) * h);
            }
            if (x < 0.0) cs = -cs; //Use antisymmetry.
            return cs;
        }
        #endregion

        public static class BesselClass
        {
            const double xj00 = 5.783185962946785,
                         xj10 = 3.047126234366209e1,
                         xj01 = 1.468197064212389e1,
                         xj11 = 4.921845632169460e1,
                         twoopi = 0.6366197723675813,
                         pio4 = 0.7853981633974483;

            static double[] j0r  = { 1.682397144220462e-4, 2.058861258868952e-5, 5.288947320067750e-7, 5.557173907680151e-9, 2.865540042042604e-11, 7.398972674152181e-14, 7.925088479679688e-17 },
                            j0s  = { 1.0, 1.019685405805929e-2, 5.130296867064666e-5, 1.659702063950243e-7, 3.728997574317067e-10, 5.709292619977798e-13, 4.932979170744996e-16 },
                            j0pn = { 9.999999999999999e-1, 1.039698629715637, 2.576910172633398e-1, 1.504152485749669e-2, 1.052598413585270e-4 },
                            j0pd = { 1.0, 1.040797262528109, 2.588070904043728e-1, 1.529954477721284e-2, 1.168931211650012e-4 },
                            j0qn = { -1.562499999999992e-2, -1.920039317065641e-2, -5.827951791963418e-3, -4.372674978482726e-4, -3.895839560412374e-6 },
                            j0qd = { 1.0, 1.237980436358390, 3.838793938147116e-1, 3.100323481550864e-2, 4.165515825072393e-4 },
                            j1r  = { 7.309637831891357e-5, 3.551248884746503e-6, 5.820673901730427e-8, 4.500650342170622e-10, 1.831596352149641e-12, 3.891583573305035e-15, 3.524978592527982e-18 },
                            j1s  = { 1.0, 9.398354768446072e-3, 4.328946737100230e-5, 1.271526296341915e-7, 2.566305357932989e-10, 3.477378203574266e-13, 2.593535427519985e-16 },
                            j1pn = { 1.0, 1.014039111045313, 2.426762348629863e-1, 1.350308200342000e-2, 9.516522033988099e-5 },
                            j1pd = { 1.0, 1.012208056357845, 2.408580305488938e-1, 1.309511056184273e-2, 7.746422941504713e-5 },
                            j1qn = { 4.687499999999991e-2, 5.652407388406023e-2, 1.676531273460512e-2, 1.231216817715814e-3, 1.178364381441801e-5 },
                            j1qd = { 1.0, 1.210119370463693, 3.626494789275638e-1, 2.761695824829316e-2, 3.240517192670181e-4 },
                            y0r  = { -7.653778457189104e-3, -5.854760129990403e-2, 3.720671300654721e-4, 3.313722284628089e-5, 4.247761237036536e-8, -4.134562661019613e-9, -3.382190331837473e-11, -1.017764126587862e-13, -1.107646382675456e-16 },
                            y0s  = { 1.0, 1.125494540257841e-2, 6.427210537081400e-5, 2.462520624294959e-7, 7.029372432344291e-10, 1.560784108184928e-12, 2.702374957564761e-15, 3.468496737915257e-18, 2.716600180811817e-21 },
                            y0pn = { 9.999999999999999e-1, 1.039698629715637, 2.576910172633398e-1, 1.504152485749669e-2, 1.052598413585270e-4 },
                            y0pd = { 1.0, 1.040797262528109, 2.588070904043728e-1, 1.529954477721284e-2, 1.168931211650012e-4 },
                            y0qn = { -1.562499999999992e-2, -1.920039317065641e-2, -5.827951791963418e-3, -4.372674978482726e-4, -3.895839560412374e-6 },
                            y0qd = { 1.0, 1.237980436358390, 3.838793938147116e-1, 3.100323481550864e-2, 4.165515825072393e-4 },
                            y1r  = { -1.041835425863234e-1, -1.135093963908952e-5, 2.212118520638132e-4, 1.270981874287763e-6, -3.982892100836748e-8, -4.820712110115943e-10, -1.929392690596969e-12, -2.725259514545605e-15 },
                            y1s  = { 1.0, 1.186694184425838e-2, 7.121205411175519e-5, 2.847142454085055e-7, 8.364240962784899e-10, 1.858128283833724e-12, 3.018846060781846e-15, 3.015798735815980e-18 },
                            y1pn = { 1.0, 1.014039111045313, 2.426762348629863e-1, 1.350308200342000e-2, 9.516522033988099e-5 },
                            y1pd = { 1.0, 1.012208056357845, 2.408580305488938e-1, 1.309511056184273e-2, 7.746422941504713e-5 },
                            y1qn = { 4.687499999999991e-2, 5.652407388406023e-2, 1.676531273460512e-2, 1.231216817715814e-3, 1.178364381441801e-5 },
                            y1qd = { 1.0, 1.210119370463693, 3.626494789275638e-1, 2.761695824829316e-2, 3.240517192670181e-4 },
                            i0p  = { 9.999999999999997e-1, 2.466405579426905e-1, 1.478980363444585e-2, 3.826993559940360e-4, 5.395676869878828e-6, 4.700912200921704e-8, 2.733894920915608e-10, 1.115830108455192e-12,
                                     3.301093025084127e-15, 7.209167098020555e-18, 1.166898488777214e-20, 1.378948246502109e-23, 1.124884061857506e-26, 5.498556929587117e-30 },
                            i0q  = { 4.463598170691436e-1, 1.702205745042606e-3, 2.792125684538934e-6, 2.369902034785866e-9, 8.965900179621208e-13 },
                            i0pp = { 1.192273748120670e-1, 1.947452015979746e-1, 7.629241821600588e-2, 8.474903580801549e-3, 2.023821945835647e-4 },
                            i0qq = { 2.962898424533095e-1, 4.866115913196384e-1, 1.938352806477617e-1, 2.261671093400046e-2, 6.450448095075585e-4, 1.529835782400450e-6 },
                            i1p  = { 5.000000000000000e-1, 6.090824836578078e-2, 2.407288574545340e-3, 4.622311145544158e-5,5.161743818147913e-7, 3.712362374847555e-9, 1.833983433811517e-11, 6.493125133990706e-14,
                                     1.693074927497696e-16, 3.299609473102338e-19, 4.813071975603122e-22, 5.164275442089090e-25, 3.846870021788629e-28, 1.712948291408736e-31 },
                            i1q  = { 4.665973211630446e-1, 1.677754477613006e-3, 2.583049634689725e-6, 2.045930934253556e-9, 7.166133240195285e-13 },
                            i1pp = { 1.286515211317124e-1, 1.930915272916783e-1, 6.965689298161343e-2, 7.345978783504595e-3, 1.963602129240502e-4 },
                            i1qq = { 3.309385098860755e-1, 4.878218424097628e-1, 1.663088501568696e-1, 1.473541892809522e-2, 1.964131438571051e-4, -1.034524660214173e-6 },
                            k0pi = { 1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6 },
                            k0qi = { 9.847324170755358e-1, 1.518396076767770e-2, 8.362215678646257e-5 },
                            k0p  = { 1.159315156584126e-1, 2.770731240515333e-1, 2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6 },
                            k0q  = { 9.836249671709183e-1, 1.627693622304549e-2, 9.809660603621949e-5 },
                            k0pp = { 1.253314137315499, 1.475731032429900e1, 6.123767403223466e1, 1.121012633939949e2, 9.285288485892228e1, 3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2 },
                            k0qq = { 1.0, 1.189963006673403e1, 5.027773590829784e1, 9.496513373427093e1, 8.318077493230258e1, 3.181399777449301e1, 4.443672926432041, 1.408295601966600e-1 },
                            k1pi = { 0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7 },
                            k1qi = { 9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5 },
                            k1p  = { -3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7 },
                            k1q  = { 9.861813171751389e-1, 1.375094061153160e-2, 6.774221332947002e-5 },
                            k1pp = { 1.253314137315502, 1.457171340220454e1, 6.063161173098803e1, 1.147386690867892e2, 1.040442011439181e2, 4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1 },
                            k1qq = { 1.0, 1.125154514806458e1, 4.427488496597630e1, 7.616113213117645e1, 5.863377227890893e1, 1.850303673841586e1, 1.857244676566022, 2.538540887654872e-2 };

            static Nmbr nump, denp, numq, denq, ans, y, z, ax, xx, j0x, j1x, logm1 = new Cmplx(0, Math.PI);

            #region private

            static void rat(Nmbr x, double[] r, double[] s, int n)
            {
                //Evaluates rational approximation.
                y = x * x;
                z = 64.0 - y;
                nump = r[n];
                denp = s[n];
                for (int i = n - 1; i >= 0; i--)
                {
                    nump = nump * z + r[i];
                    denp = denp * y + s[i];
                }
            }

            static void asp(double[] pn, double[] pd, double[] qn, double[] qd, double fac)
            {
                //Evaluates asymptotic approximation.
                z = 8.0 / ax;
                y = z * z;
                xx = ax - fac * pio4;
                nump = pn[4];
                denp = pd[4];
                numq = qn[4];
                denq = qd[4];
                for (int i = 3; i >= 0; i--)
                {
                    nump = nump * y + pn[i];
                    denp = denp * y + pd[i];
                    numq = numq * y + qn[i];
                    denq = denq * y + qd[i];
                }
            }

            static Nmbr poly(double[] cof, int n, Nmbr x)
            {
                //Evaluate a polynomial.
                ans = cof[n];
                for (int i = n - 1; i >= 0; i--) ans = ans * x + cof[i];
                return ans;
            }

            static Nmbr j0(Nmbr x)
            {
                //Returns the Bessel functionJ0.x/for any real x.
                if (Maths.Abs(x) < 8.0)
                {
                    //Direct rational function fit.
                    rat(x, j0r, j0s, 6);
                    return nump * (y - xj00) * (y - xj10) / denp;
                }
                else
                {
                    //Fitting function (6.5.9).
                    ax = x.x < 0 ? -x : x; asp(j0pn, j0pd, j0qn, j0qd, 1.0);
                    return Maths.Sqrt(twoopi / ax) * (Maths.Cos(xx) * nump / denp - z * Maths.Sin(xx) * numq / denq);
                }
            }

            static Nmbr j1(Nmbr x)
            {
                //Returns the Bessel functionJ1.x/for any real x.
                if (Maths.Abs(x) < 8.0)
                {
                    //Direct rational approximation.
                    rat(x, j1r, j1s, 6);
                    return x * nump * (y - xj01) * (y - xj11) / denp;
                }
                else
                {
                    //Fitting function (6.5.9).
                    ax = x.x < 0 ? -x : x; asp(j1pn, j1pd, j1qn, j1qd, 3.0);
                    Nmbr ans = Maths.Sqrt(twoopi / ax) * (Maths.Cos(xx) * nump / denp - z * Maths.Sin(xx) * numq / denq);
                    return x.x >= 0.0 ? ans : -ans;
                }
            }

            static Nmbr y0(Nmbr x)
            {
                //Returns the Bessel function Y0(x)/for positive x.
                j0x = j0(x);
                if (Maths.Abs(x) < 8.0)
                {
                    //Rational function approximation of (6.5.8).
                    rat(x, y0r, y0s, 8);
                    return nump / denp + twoopi * j0x * Maths.Log(x);
                }
                else
                {
                    //Fitting function (6.5.10).
                    asp(y0pn, y0pd, y0qn, y0qd, 1.0);
                    ans = Maths.Sqrt(twoopi / ax) * (Maths.Sin(xx) * nump / denp + z * Maths.Cos(xx) * numq / denq);
                    if(x.x < 0 && x.y == 0) ans += twoopi * j0x * logm1;
                    return ans;
                }
            }

            static Nmbr y1(Nmbr x)
            {
                //Returns the Bessel functionY1.x/for positive x.
                j1x = j1(x);
                if (Maths.Abs(x) < 8.0)
                {
                    //Rational function approximation of (6.5.8).
                    rat(x, y1r, y1s, 7);
                    return x * nump / denp + twoopi * (j1x * Maths.Log(x) - 1.0 / x);
                }
                else
                {
                    //Fitting function (6.5.10).
                    asp(y1pn, y1pd, y1qn, y1qd, 3.0);
                    ans = Maths.Sqrt(twoopi / ax) * (Maths.Sin(xx) * nump / denp + z * Maths.Cos(xx) * numq / denq);
                    if (x.x < 0 && x.y == 0) ans += twoopi * j1x * logm1;
                    return ans;
                }
            }

            static Nmbr i0(Nmbr x)
            {
                //Returns the Bessel function I0(x)/for any real x.
                if (Maths.Abs(x) < 15.0)
                {
                    //Rational approximation.
                    y = x * x;
                    return poly(i0p, 13, y) / poly(i0q, 4, 225.0 - y);
                }
                else
                {
                    //Rational approximation with
                    ax = x.x < 0 ? -x : x; z = 1.0 - 15.0 / ax;
                    return Maths.Exp(ax) * poly(i0pp, 4, z) / (poly(i0qq, 5, z) * Maths.Sqrt(ax));
                }
            }

            static Nmbr i1(Nmbr x)
            {
                //Returns the Bessel functionJ1(x)for any real x.
                if ((ax = Maths.Abs(x.x)) < 15.0)
                {
                    //Rational approximation.
                    y = x * x;
                    return x * poly(i1p, 13, y) / poly(i1q, 4, 225.0 - y);
                }
                else
                {
                    //Rational approximation with
                    ax = x.x < 0 ? -x : x; z = 1.0 - 15.0 / ax;
                    Nmbr ans = Maths.Exp(ax) * poly(i1pp, 4, z) / (poly(i1qq, 5, z) * Maths.Sqrt(ax));
                    return x.x > 0.0 ? ans : -ans;
                }
            }

            static Nmbr k0(Nmbr x)
            {
                //Returns the Bessel functionY0.x/for positive x.
                if (x <= 1.0)
                {
                    //Rational function approximation of (6.5.8).
                    z = x * x;
                    Nmbr term = poly(k0pi, 4, z) * Maths.Log(x) / poly(k0qi, 2, 1.0 - z);
                    return poly(k0p, 4, z) / poly(k0q, 2, 1.0 - z) - term;
                }
                else
                {
                    // Rational approxim
                    z = 1.0 / x;
                    return Maths.Exp(-x) * poly(k0pp, 7, z) / (poly(k0qq, 7, z) * Maths.Sqrt(x));
                }
            }

            static Nmbr k1(Nmbr x)
            {
                //Returns the Bessel function K1(x)for positive x.
                if (x <= 1.0)
                {
                    //Rational function approximation of (6.5.8).
                    z = x * x;
                    Nmbr term = poly(k1pi, 4, z) * Maths.Log(x) / poly(k1qi, 2, 1.0 - z);
                    return x * (poly(k1p, 4, z) / poly(k1q, 2, 1.0 - z) + term) + 1.0 / x;
                }
                else
                {
                    //Fitting function (6.5.10).
                    z = 1.0 / x;
                    return Maths.Exp(-x) * poly(k1pp, 7, z) / (poly(k1qq, 7, z) * Maths.Sqrt(x));
                }
            }

            #endregion

            public static Nmbr Jn(int n, Nmbr x)
            {
                Nmbr ACC = 100, ans, sum, bj, bjm, bjp, tox, BIGNO = 1.0e10, BIGNI = 1.0e-10; bool jsum; int j, m;
                if (n == 0) return j0(x);
                if (n == 1) return j1(x);
                ax = Maths.Abs(x);
                if (ax * ax < 8 * 1e-15) return 0;
                else if (ax > n)
                {
                    tox = 2.0 / x;
                    bj = j1(x); //Starting values for the recurrence.
                    bjm = j0(x);
                    for (j = 1; j < n; j++)
                    {  //Recurrence (6.5.7).
                        bjp = j * tox * bj - bjm;
                        bjm = bj;
                        bj = bjp;
                    }
                    ans = bj;
                }
                else
                {
                    tox = 2.0 / x;
                    m = 2 * ((n + (int)Maths.Sqrt(ACC * n)) / 2);
                    jsum = false;
                    bjp = ans = sum = 0.0;
                    bj = 1.0;
                    for (j = m; j > 0; j--)
                    {
                        //The downward recurrence.
                        bjm = j * tox * bj - bjp;
                        bjp = bj;
                        bj = bjm;
                        if (Maths.Abs(bj) > BIGNO)
                        { //Renormalize to prevent overflows.
                            bj *= BIGNI;
                            bjp *= BIGNI;
                            ans *= BIGNI;
                            sum *= BIGNI;
                        }
                        if (jsum) sum += bj; //Accumulate the sum.
                        jsum = !jsum; //Change false to true or vice versa.
                        if (j == n) ans = bjp; //Save the unnormalized answer.
                    }
                    sum = 2.0 * sum - bj; //Compute (5.4.16)
                    ans /= sum;     //and use it to normalize the answer.
                }
                return ans;
            }

            public static Nmbr Yn(int n, Nmbr x)
            {
                Nmbr by, bym, byp, tox; int j;
                if (n == 0) return y0(x);
                if (n == 1) return y1(x);
                tox = 2.0 / x;
                by = y1(x); //Starting values for the recurrence.
                bym = y0(x);
                for (j = 1; j < n; j++)
                {  //Recurrence (6.5.7).
                    byp = j * tox * by - bym;
                    bym = by;
                    by = byp;
                }
                return by;
            }

            public static Nmbr In(int n, Nmbr x)
            {
                Nmbr ACC = 100, BIGNO = 1.0e10, BIGNI = 1.0e-10, ans, bi, bim, bip, tox; int j, m;
                if (n == 0) return i0(x);
                if (n == 1) return i1(x);
                ax = Maths.Abs(x);
                if (ax * ax < 8 * double.MinValue) return 0;
                else
                {
                    tox = 2.0 / x;
                    m = 2 * ((n + (int)(Maths.Sqrt(ACC * n))) / 2);
                    bip = ans = 0.0;
                    bi = 1.0;
                    for (j = m; j > 0; j--)
                    {
                        //The downward recurrence.
                        bim = j * tox * bi + bip;
                        bip = bi;
                        bi = bim;
                        if (Maths.Abs(bi) > BIGNO)
                        { //Renormalize to prevent overflows.
                            bi *= BIGNI;
                            bip *= BIGNI;
                            ans *= BIGNI;
                        }
                        if (j == n) ans = bip; //Save the unnormalized answer.
                    }
                    ans *= i0(x) / bi;     //and use it to normalize the answer.
                }
                return ans;
            }

            public static Nmbr Kn(int n, Nmbr x)
            {
                Nmbr bk, bkm, bkp, tox; int j;
                if (n == 0) return k0(x);
                if (n == 1) return k1(x);
                tox = 2.0 / x;
                bk = k1(x); //Starting values for the recurrence.
                bkm = k0(x);
                for (j = 1; j < n; j++)
                {  //Recurrence (6.5.7).
                    bkp = j * tox * bk + bkm;
                    bkm = bk;
                    bk = bkp;
                }
                return bk;
            }
        }
    }

    class SpecialFunctionsException : Exception
    {
        /// <summary>
        /// Exception from the Matrix Class
        /// </summary>
        /// <param name="Message">Message to be displayed</param>
        public SpecialFunctionsException(string Message)
            : base(Message)
        { }
    }
}
