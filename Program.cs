using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;

namespace WissRech
{
    internal class Program
    {
        static void Main(string[] args)
        {
            //       "Hilbert"

            //"Identity"

            // "FE"


            int N = 3;
            double[] b = new double[N];
            double[] txk = new double[N];
            for (int i = 0; i < N; i++)
            {
                b[i] = 1;
                txk[i] = 1;
            }

            
            double[] x = {0,0,0,1,1,1,0,0,1,1,0,1,1,0,1,1,2,1,1,0,2,0,2,1,0,1,1,1,1,2,0,1,1,2,0,2,1,1,2,1,2,2,1,1,2,2,1,2};
            double[] N_B = { 0, 1, 0, 2, 2, 1, 2, 2, 1, 2};
            double[] D_B = { 0, 1 ,0,0,0,2};
            double[] z = { 1,2,1,2,1,2,1,2,1};
            FEM D = new FEM(x, 2,N_B,D_B);
            D.print_points();
            D.print_global_indeces();
            D.print_set_i_N();
            D.print_set_D();
            

            for(int i = 0; i < 9; i++)
            {
                //Console.WriteLine(x.Length/6);
                Console.WriteLine("ERGEBNIS" + D.stiffness_vec(z)[i]);
            }
            //D.stiffness_vec(z);
            Console.ReadKey();
        }
    
    }
}
