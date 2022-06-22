using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SciSharp
{
    public class Indexer
    {
        int[] indx;
        public Indexer() { }

        public Indexer(int start, int end)
        {
            List<int> i = new List<int> { start++ };
            while (start <= end) i.Add(start++);
            indx = i.ToArray();
        }

        public Indexer(int start, int step, int end)
        {
            List<int> i = new List<int> { start };
            start += step;
            while (start <= end)
            { i.Add(start); start += step; }
            indx = i.ToArray();
        }

        public Indexer(int[] array)
        { indx = array.ToArray(); }

        public static implicit operator int[](Indexer c) => c.indx;
        public static implicit operator Indexer(int[] c) => new Indexer(c);

        // OPERATORS
        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator +(Indexer c1, Indexer c2)
        {
            if (c1.indx.Length != c2.indx.Length)
                throw new Exception("Indexers must be of the same length");
            Indexer c = new Indexer(); c.indx = new int[c1.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1.indx[i] + c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator +(int c1, Indexer c2)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 + c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator +(Indexer c2, int c1)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 + c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator -(Indexer c1, Indexer c2)
        {
            if (c1.indx.Length != c2.indx.Length)
                throw new Exception("Indexers must be of the same length");
            Indexer c = new Indexer(); c.indx = new int[c1.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1.indx[i] - c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator -(int c1, Indexer c2)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 - c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator -(Indexer c2, int c1)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c2.indx[i] - c1; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Indexer operator -(Indexer c)
        {
            Indexer b = new Indexer(); c.indx = new int[c.indx.Length];
            for (int i = 0; i < b.indx.Length; i++)
            { b.indx[i] = -c.indx[i]; }
            return b;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator *(Indexer c1, Indexer c2)
        {
            if (c1.indx.Length != c2.indx.Length)
                throw new Exception("Indexers must be of the same length");
            Indexer c = new Indexer(); c.indx = new int[c1.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1.indx[i] * c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator *(int c1, Indexer c2)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 * c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator *(Indexer c2, int c1)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 * c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator /(Indexer c1, Indexer c2)
        {
            if (c1.indx.Length != c2.indx.Length)
                throw new Exception("Indexers must be of the same length");
            Indexer c = new Indexer(); c.indx = new int[c1.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1.indx[i] / c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator /(int c1, Indexer c2)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c1 / c2.indx[i]; }
            return c;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Indexer operator /(Indexer c2, int c1)
        {
            Indexer c = new Indexer(); c.indx = new int[c2.indx.Length];
            for (int i = 0; i < c.indx.Length; i++)
            { c.indx[i] = c2.indx[i]/c1; }
            return c;
        }
    }
}
