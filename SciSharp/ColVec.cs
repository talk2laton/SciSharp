using System;
using System.Collections.Generic;
using System.Text;

namespace SciSharp
{
    public class ColVec:Matrix
    {
        public ColVec(Nmbr[] entry)
        {
            rows = entry.Length;
            cols = 1;
            mat = new List<List<Nmbr>>(cols);
            for (int i = 0; i < rows; i++)
            {
                mat.Add(new List<Nmbr>(cols));
                for (int j = 0; j < cols; j++)
                {
                    mat[i][j] = entry[i];
                }
            }
        }


    }
}
