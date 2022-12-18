using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM_final
{
    internal class Function
    {
        public Function()
        {

        }

        public double get_f(double[]x)
        {
            double outp;
            outp = x[1] + x[0];
            return outp;

        }

     
    }
}

