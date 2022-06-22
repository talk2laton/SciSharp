using System;
using System.Collections.Generic;
using System.Text;

namespace SciSharp
{
    public class RowVec:Matrix
    {
        public RowVec(Nmbr[] entry)
        {
            rows = 1;
            cols = entry.Length;
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
