using System;
using SciSharp;
using System.Collections.Generic;
using System.Linq;

namespace SciSharpTester
{
    class Program
    {
        static void Main(string[] args)
        {
            Maths.format = Maths.Format.Long;
            Random rand = new Random();
            List<int> I = new List<int>();
            List<Nmbr> C = new List<Nmbr>();
            Console.Write("clc; \n a = [");
            for (int i = 0; i < 50; i++)
            {
                double x = rand.NextDouble(), y = rand.NextDouble();
                C.Add(12 * new Cmplx(x, y) + new Cmplx(5, 5));
                I.Add((int)(5 * rand.NextDouble()));
                Console.WriteLine("besseli({0}, {1})", I.Last(), C.Last());
            }
            Console.WriteLine("];");

            Console.Write("b = [");
            for (int i = 0; i < I.Count; i++)
            {
                Console.WriteLine(Maths.BesselI(I[i], C[i]));
            }
            Console.WriteLine("];");

            Console.WriteLine("disp(abs((a-b)./a));");

            //Console.WriteLine(Maths.BesselI(1, new Cmplx(3,3)));
            //Console.WriteLine(Maths.BesselI(1, new Cmplx(18, 18)));
        }
    }
}
