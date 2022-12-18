using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace FEM_final
{
    public class Vektoren
    {
        private double[] x;
        private int dimension;
        public Vektoren(double[] x)
        {
            dimension = x.Length;
            objekt();
            this.x = x;
        }
        public double[] get_vektor()
        {
            return x;
        }
        private void objekt()
        {
            x = new double[dimension];
        }
     
        public void AddElement(double o)
        {
            Array.Resize(ref x, dimension + 1);
            x[dimension] = o;
            dimension = dimension + 1;


        }
        public void multiply(double c)
        {
            for(int j = 0; j < dimension; j++)
            {
                x[j]=c*x[j];
            }
        }
        public void Ausgabe()
        {
            for(int j = 0; j < dimension; j++)
            {
                Console.WriteLine(x[j]);
            }
        }
        public Vektoren add(Vektoren y)
        {
           
            if (x.Length!=y.get_dimension())
            {
                throw new InvalidOperationException("Dimensionen passen nicht");
            }
            else
            {
                double[] yi = new double[dimension];
                yi = y.get_vektor();
                double[] outp = new double[dimension];

                for (int j = 0;j<dimension;j++)
                {
                    outp[j] = x[j] + yi[j];
                }
                Vektoren output = new Vektoren(outp);
                return output;
            }
            
        }
        public int get_dimension()
        {
            return dimension;
        }
    }
}
