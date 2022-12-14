using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Net.NetworkInformation;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;

namespace WissRech
{
    public class FEM
    {
        private double[] T;
        private int N;
        private double [] N_boundary;
        private double[] D_boundary;
        private double[] points; // globale Knotenpunkte (x_1,x_2)
        private int[] global_index; // globaler Index zu jedem Knotenpunkt eines Dreiecks
        private int[] N_i_index;
        private int[] D_i_index;
        public int K;
        private int[] r_Wert;
        private int[] k_Wert;
 

        /// <summary>
        /// 
        /// </summary>
        /// <param name="T">Triangulierung</param>
        /// <param name="N">Dimension</param>
        /// <param name="N_boundary">Neumann Boundary</param>
        /// <param name="D_boundary">Dirichlet boundary</param>
        public FEM(double[] T, int N, double[] N_boundary, double[] D_boundary) //N={1,2}, T speichert die Punkte der Triangulation, 1D: 2N, KONSTRUKTOR
        {
            this.T = T;
            this.N = N;
            this.N_boundary=N_boundary;
            this.D_boundary=D_boundary;
            K = T.Length / 6;
            points = point(T);
            global_index = global_indices(T);
            N_i_index = set_i_N(T);
            D_i_index = set_D(global_index,N_i_index);
            r_Wert = r_Werte(T);
            k_Wert = k_Werte(T);
        }

        public double[] point(double[] T) // legt Vektor mit Knoten der Triangulation an
        {
            double[] points_1 = new double[2];
            points_1[0] = T[0];
            points_1[1] = T[1];
            for (int i = 2; i < T.Length; i = i + 2)
            {
                int k = 0;
                for (int j = 0; j < points_1.Length; j = j + 2)
                {
                    if (T[i] == points_1[j] & T[i + 1] == points_1[j + 1])
                    {
                        k++;
                    }
                }
                if (k == 0)
                {
                    Array.Resize(ref points_1, points_1.Length + 2);
                    points_1[points_1.Length - 2] = T[i];
                    points_1[points_1.Length - 1] = T[i + 1];
                }
            }
            return points_1;
        }

        public int[] global_indices(double[] T)
        {
            int[] global_indices = new int[T.Length/2];
            for( int i=0; i < T.Length; i += 2)
            {
                for( int j=0; j < points.Length; j+=2)
                {
                    if (T[i] == points[j] & T[i + 1] == points[j + 1])
                        {
                        global_indices[i/2] = j/2;
                        }
                }
            }
            return(global_indices);
        }

        public int[] set_i_N(double[] T)
        {
            int[] set_i_N = new int[0];
            for( int i=0; i<T.Length; i += 2)
            {
                int k = 0;
                for (int j = 0; j < D_boundary.Length; j += 2)
                {
                    if (T[i] == D_boundary[j] & T[i + 1] == D_boundary[j + 1])
                    {
                        k++;
                    }
                }
                 
                if(k==0) 
                {
                    int l = 0;
                    for (int m = 0; m < set_i_N.Length; m++)
                    {
                        if (global_index[i / 2] == set_i_N[m])
                        {
                            l++;
                        }
                    }
                    if (l == 0)
                    {
                        Array.Resize(ref set_i_N, set_i_N.Length + 1);
                        set_i_N[set_i_N.Length - 1] = global_index[i / 2];
                    }
                }

            }
            return set_i_N;
        }

        public int [] set_D(int[] global_index, int[] N_i_index)
        {
            int[] set_D = new int[0];
            for(int i=0; i<global_index.Length; i++)
            {
                int m = 0;
                for (int j = 0; j < N_i_index.Length; j++)
                {
                    if (global_index[i] == N_i_index[j])
                    {
                        m++;
                    }
                    for (int k = 0; k < set_D.Length; k++)
                    {
                        if (set_D[k] == global_index[i])
                        {
                            m++;
                        }
                    }
                }
                if (m == 0)
                {
                    Array.Resize(ref set_D, set_D.Length + 1);
                    set_D[set_D.Length - 1] = global_index[i];
                }
            }
            return set_D;
        }

       public double A_r_s(int k,int r, int s)
        {
            double det_F_k = (T[6*k + 2] - T[6*k]) * (T[6*k + 5] - T[6*k + 1]) - (T[6*k + 3] - T[6*k + 1]) * (T[6*k + 4] - T[6*k]);
            double erg_integral = 0.5 * Math.Abs(det_F_k); // 1/2 ist Flächeninhalt von T_0, |det_F_k| kommt aus Trafo
            double[] erg1 = new double[2]; //r=0, drei Gradienten aus dem Integrant
            double[] erg2 = new double[2];  //r=1
            double[] erg3 = new double[2]; //r=2
            erg1[0] = 1 / det_F_k * (T[k + 3] - T[k + 5]);
            erg1[1] = 1 / det_F_k * (T[k + 4] - T[k+2]);
            erg2[0] = 1 / det_F_k * (T[k + 5] - T[k + 1]);
            erg2[1] = 1 / det_F_k * (T[k] - T[k + 4]);
            erg3[0] = 1 / det_F_k * (T[k + 1] - T[k + 3]);
            erg3[1] = 1 / det_F_k * (T[k + 2] - T[k]);
            //Console.WriteLine("erg[i]"+erg1[1]);
            //Console.WriteLine("erg_integral"+ erg_integral);
            //Console.WriteLine("ABS" + Math.Abs(det_F_k));
            //Console.ReadKey();

            if(r==0 & s == 0)
            {
                erg_integral = (erg1[0] * erg1[0] + erg1[1] * erg1[1]) * erg_integral;
                //Console.WriteLine("Yes");
                //Console.WriteLine(erg_integral);
            }
            if(r==0 & s==1 | s==0 & r == 1)
            {
                erg_integral *= (erg1[0] * erg2[0] + erg1[1] * erg2[1]);
            }
            if (r == 0 & s == 2 | r==2 & s == 0)
            {
                erg_integral *= (erg1[0] * erg3[0] + erg1[1] * erg3[1]);
            }
            if(r==1 & s == 1)
            {
                erg_integral *= (erg2[0] * erg2[0] + erg2[1] * erg2[1]);
            }
            if(r==1 & s==2 | s==1 & r == 2)
            {
                erg_integral *= (erg2[0] * erg3[0] + erg2[1] * erg3[1]);
            }
            if(r==2 & s == 2)
            {
                erg_integral*= (erg3[0] * erg3[0] + erg3[1] * erg3[1]);
            }

            return erg_integral;

        }

        public double[] stiffness_vec(double[] x)
        {
            double[] erg= new double[x.Length]; //Ergebnisvektor von Matrix-Vektor-Multiplikation, Länge x= Anzahl Knoten
            //Console.WriteLine("Länge x: "+x.Length);
            //Console.ReadKey();
            for( int i=0; i < x.Length; i++) // (Ax)_i=erg_i
            {
                for(int j=0; j < x.Length; j++)
                {
                    double A_ij = 0; //Eintrag A_ij in Matrix
                    for( int k=0; k<K; k++) //läuft über Dreiecke
                    {
                        for(int r = 0; r < 3; r++)
                        {
                            for (int s = 0; s < 3; s++)
                            {
                                if(global_index[(6*k+2*r)/2]==i & global_index[(6*k+2*s)/2]==j) //nach der for-Schleife ist Eintrag A_ij fertig
                                {
                                    A_ij+= A_r_s(k, r, s);
                                }
                            }
                        }
                    }
                    erg[i] += A_ij * x[j];
                }
            }
                return erg;
        }

        public int get_r(int k) //liefert zu einer Numerierung k von einem Dreieck ein r-Wert, st. i(r,h)=k für irgendein h.
         { 
            for(int i =0; i < T.Length/2; i++)
            {
                if (global_index[i] == k)
                {
                    return r_Wert[i];
                     
                }
            }
            return 0;
        }

        public int[] r_Werte(double[] T)
        {
            int[] r_Werte = new int[T.Length / 2];
            for(int i=0; i<r_Werte.Length; i+=3)
            {
                r_Werte[i] = 0;
                r_Werte[i + 1] = 1;
                r_Werte[i + 2] = 2;
            }
            return r_Werte; 
        }

        public int [] k_Werte(double[] T)
        {
            int[] k_Werte = new int[T.Length / 6];
            for(int i=0; i< k_Werte.Length; i++)
            {
                k_Werte[i] = i;
            }
            return k_Werte;
        }

 

        public void print_points()
        {
            for(int i=0;i<points.Length; i+=2)
            {
                Console.WriteLine("("+points[i]+"," +points[i+1]+")");

            }
        }
        public void print_global_indeces()
        {
            for (int i = 0; i < global_index.Length; i++)
            {
                Console.WriteLine("Punkt " + "("+T[2 * i]+","+T[2*i+1]+") hat Index " + global_index[i] + " ,r_Wert: " + r_Wert[i]+"" +
                    ", k_Wert: " + k_Wert[2*i/6]);

            }

        }

        public void print_set_i_N()
        {
            Console.WriteLine("Neumann + Interior");
            for (int i = 0; i < N_i_index.Length; i++)
            {
                Console.WriteLine("Index "+N_i_index[i]);
            }

        }
        public void print_set_D()
        {
            Console.WriteLine("Dirichlet");
            for (int i = 0; i < D_i_index.Length; i++)
            {
                Console.WriteLine("Index " + D_i_index[i]);
            }

        }

        public void number_triangles()
        {
            Console.WriteLine(K);
        }

    }

}
