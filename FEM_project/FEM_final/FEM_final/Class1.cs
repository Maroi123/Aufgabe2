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

namespace FEM_final
{
    public class FEM
    {
        private double[] T;
        private int N;
        private double[] N_boundary;
        private double[] D_boundary;
        private double[] points; // globale Knotenpunkte (x_1,x_2)
        private int[] global_index; // globaler Index zu jedem Knotenpunkt eines Dreiecks
        private int[] N_i_index;
        private int[] D_i_index;
        public int K;
        private int[] r_Wert;
        private int[] k_Wert;
        Function f;



        public FEM(double[] T, int N, double[] N_boundary, double[] D_boundary) //N={1,2}, T speichert die Punkte der Triangulation, 1D: 2N, KONSTRUKTOR
        {
            this.T = T;
            this.N = N;
            this.N_boundary = N_boundary;
            this.D_boundary = D_boundary;
            K = T.Length / 6;
            points = point(T);
            global_index = global_indices(T);
            N_i_index = set_i_N(T);
            D_i_index = set_D(global_index, N_i_index);
            r_Wert = r_Werte(T);
            k_Wert = k_Werte(T);
            f = new Function();
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
            int[] global_indices = new int[T.Length / 2];
            for (int i = 0; i < T.Length; i += 2)
            {
                for (int j = 0; j < points.Length; j += 2)
                {
                    if (T[i] == points[j] & T[i + 1] == points[j + 1])
                    {
                        global_indices[i / 2] = j / 2;
                    }
                }
            }
            return (global_indices);
        }

        public int[] set_i_N(double[] T)
        {
            int[] set_i_N = new int[0];
            for (int i = 0; i < T.Length; i += 2)
            {
                int k = 0;
                for (int j = 0; j < D_boundary.Length; j += 2)
                {
                    if (T[i] == D_boundary[j] & T[i + 1] == D_boundary[j + 1])
                    {
                        k++;
                    }
                }

                if (k == 0)
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

        public int[] set_D(int[] global_index, int[] N_i_index)
        {
            int[] set_D = new int[0];
            for (int i = 0; i < global_index.Length; i++)
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

        public double A_r_s(int k, int r, int s)
        {
            double det_F_k = (T[6 * k + 2] - T[6 * k]) * (T[6 * k + 5] - T[6 * k + 1]) - (T[6 * k + 3] - T[6 * k + 1]) * (T[6 * k + 4] - T[6 * k]);
            double erg_integral = 0.5 * Math.Abs(det_F_k); // 1/2 ist Flächeninhalt von T_0, |det_F_k| kommt aus Trafo
            double[] erg1 = new double[2]; //r=0, drei Gradienten aus dem Integrant
            double[] erg2 = new double[2];  //r=1
            double[] erg3 = new double[2]; //r=2
            erg1[0] = 1 / det_F_k * (T[k + 3] - T[k + 5]);
            erg1[1] = 1 / det_F_k * (T[k + 4] - T[k + 2]);
            erg2[0] = 1 / det_F_k * (T[k + 5] - T[k + 1]);
            erg2[1] = 1 / det_F_k * (T[k] - T[k + 4]);
            erg3[0] = 1 / det_F_k * (T[k + 1] - T[k + 3]);
            erg3[1] = 1 / det_F_k * (T[k + 2] - T[k]);

            if (r == 0 & s == 0)
            {
                erg_integral = (erg1[0] * erg1[0] + erg1[1] * erg1[1]) * erg_integral;
            }
            if (r == 0 & s == 1 | s == 0 & r == 1)
            {
                erg_integral *= (erg1[0] * erg2[0] + erg1[1] * erg2[1]);
            }
            if (r == 0 & s == 2 | r == 2 & s == 0)
            {
                erg_integral *= (erg1[0] * erg3[0] + erg1[1] * erg3[1]);
            }
            if (r == 1 & s == 1)
            {
                erg_integral *= (erg2[0] * erg2[0] + erg2[1] * erg2[1]);
            }
            if (r == 1 & s == 2 | s == 1 & r == 2)
            {
                erg_integral *= (erg2[0] * erg3[0] + erg2[1] * erg3[1]);
            }
            if (r == 2 & s == 2)
            {
                erg_integral *= (erg3[0] * erg3[0] + erg3[1] * erg3[1]);
            }

            return erg_integral;

        }
        /// <summary>
        /// Rechnet das Produkt von vektoren aus
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        private double produkt(double[] x, double[] y)
        {
            double outp = 0;
            if (x.Length != y.Length)
            {
                throw new InvalidOperationException("Dimensionen passen nicht");

            }
            else
            {
                for (int i = 0; i < x.Length; i++)
                {

                    outp = outp + x[i] * y[i];
                }
                return outp;
            }
        }
        /// <summary>
        /// Implementiert das CG verfahren für die Matrix A
        /// </summary>
        /// <param name="eps">gewünschter Fehler</param>
        /// <param name="xk">Startwert</param>
        /// <param name="f">Rechte Seite f</param>
        /// <param name="Iteration">Anzahl der maximalen iterationen</param>
        public double[] CG_method(double eps, double[] xk, double[] f, int Iteration)
        {
            double fehler;
            int dimension = xk.Length;
            int i = 0;
            double bk;
            double ak;
            double[] trk1 = new double[dimension];
            double[] tdk1 = new double[dimension];
            double[] xk1 = new double[dimension];
            double[] vector = new double[dimension];
            double[] trk = new double[dimension];
            double[] tdk = new double[dimension];

            vector = stiffness_vec(xk);
            for (int j = 0; j < dimension; j++)
            {
                trk[j] = f[j] - vector[j];
                tdk[j] = trk[j];
            }
            do
            {
                vector = stiffness_vec(tdk);
                if (produkt(tdk, vector) == 0)
                {
                    ak = 0;
                }
                else
                {
                    ak = produkt(trk, trk) / produkt(tdk, vector);
                }

                for (int j = 0; j < dimension; j++)
                {
                    xk1[j] = xk[j] + ak * tdk[j];
                    trk1[j] = trk[j] - ak * vector[j];
                }
                if (produkt(trk, trk) == 0)
                {
                    bk = 0;
                }
                else
                {
                    bk = produkt(trk1, trk1) / produkt(trk, trk);
                }

                for (int j = 0; j < dimension; j++)
                {
                    tdk1[j] = trk1[j] + bk * tdk[j];
                }

                for (int j = 0; j < dimension; j++)
                {
                    xk[j] = xk1[j];
                    tdk[j] = tdk1[j];
                    trk[j] = trk1[j];

                }
                //berechne den fehler:
                trk1 = stiffness_vec(xk);
                for (int j = 0; j < dimension; j++)
                {
                    vector[j] = f[j] - trk1[j];
                }
                fehler = produkt(vector, vector);

                i++;
            } while (i < Iteration && fehler > eps);

            return xk;


        }
        public double[] stiffness_vec(double[] x)
        {
            double[] erg = new double[x.Length]; //Ergebnisvektor von Matrix-Vektor-Multiplikation, Länge x= Anzahl Knoten
            for (int i = 0; i < x.Length; i++) // (Ax)_i=erg_i
            {
                for (int j = 0; j < x.Length; j++)
                {
                    double A_ij = 0; //Eintrag A_ij in Matrix
                    for (int k = 0; k < K; k++) //läuft über Dreiecke
                    {
                        for (int r = 0; r < 3; r++)
                        {
                            for (int s = 0; s < 3; s++)
                            {
                                if (global_index[(6 * k + 2 * r) / 2] == i & global_index[(6 * k + 2 * s) / 2] == j) //nach der for-Schleife ist Eintrag A_ij fertig
                                {
                                    A_ij += A_r_s(k, r, s);
                                }
                            }
                        }
                    }
                    erg[i] += A_ij * x[j];
                }
            }
            return erg;
        }

        public double[] local_load_vec(int k)
        {
            double[] erg = new double[3];
            double det_F_k = (T[6 * k + 2] - T[6 * k]) * (T[6 * k + 5] - T[6 * k + 1]) - (T[6 * k + 3] - T[6 * k + 1]) * (T[6 * k + 4] - T[6 * k]);
            det_F_k = Math.Abs(det_F_k);
            double x = (double)1 / (double)6;
            double[] x_1 = { x, x };
            double[] x_2 = { 4 * x, x };
            double[] x_3 = { x, 4 * x };
            erg[0] = x * det_F_k * (f.get_f(F_k(k, x_1)) * (1 - x - x) + f.get_f(F_k(k, x_2)) * (1 - 4 * x - x) + f.get_f(F_k(k, x_3)) * (1 - x - 4 * x));
            erg[1] = x * det_F_k * (f.get_f(F_k(k, x_1)) * (1 - x) + f.get_f(F_k(k, x_2)) * (1 - 4 * x) + f.get_f(F_k(k, x_3)) * (1 - x));
            return erg;
        }

        public double[] global_load_vec()
        {
            double[]erg= new double[points.Length];
            for(int j = 0; j < erg.Length; j++)
            {
                for(int k = 0; k < K; k++)
                {
                    for(int r = 0; r < 3; r++)
                    {
                        if (global_index[(6 * k + 2 * r)/2] == j)
                        {
                            erg[j] += local_load_vec(k)[r];
                        }
                    }
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

        public double[] F_k(int k, double[] x)
        {
            double[] erg = new double[2];
            erg[0] = x[0] * (T[6 * k + 2] - T[6 * k]) + x[1] * (T[6 * k + 4] - T[k]) + T[6 * k];
            erg[0] = x[0] * (T[6 * k + 3] - T[6 * k+1]) + x[1] * (T[6 * k + 5] - T[k+1]) + T[6 * k+1];
            return erg;
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
            Console.WriteLine("globale Neumann + Interior Indizes: ");
            for (int i = 0; i < N_i_index.Length; i++)
            {
                Console.WriteLine("Index "+N_i_index[i]);
            }

        }
        public void print_set_D()
        {
            Console.WriteLine("globale Dirichlet Indizes: ");
            for (int i = 0; i < D_i_index.Length; i++)
            {
                Console.WriteLine("Index " + D_i_index[i]);
            }

        }

        public void number_triangles()
        {
            Console.WriteLine("Anzahl der Dreiecke: "+K);
        }




    }

}
