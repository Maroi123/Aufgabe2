using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    internal class FEM_3012
    {
        private int[,] A;
        private double[,] points;
        private int[] D_index;
        private int M;
        private int K;
        private int N;
        private Func<double[], double> f;
        private Func<double[], double> g;
        public double[,] stiffness;
        public FEM_3012(int[,] A,double[,] points,int[] D_index, Func<double[], double>f, Func<double[], double>g) //Konstruktor
        {
            this.A = A; //enthält die globalen Indices der Dreiecke
            this.points = points; //enthält die Koordinaten der nodes
            this.D_index = D_index; //Dirichlet Rand
            M=A.GetLength(0); //Anzahl der Eckpunkte der Dreiecke, hier 3
            K=A.GetLength(1); // ANzahl der Dreiecke
            N = points.GetLength(1); //Anzahl der nodes
            this.f = f; //rechte Seite von -Laplace=f
            this.g = g; // g setzt Dirichlet boundary
            stiffness = global_stiffness();
        }
        
        public double p(int k, int r, int coord) //gibt Eckpunkt von Dreieck k, Ecke r zurück
        {
            double erg= points[coord, A[r, k]];
            return erg;
        }

        public double[,] local_stiffness(int k) //berechnet locale stiffness Matrix 
        {
            double[,] stiffness = new double[M, M]; //3x3 Matrix
            double det= 0.5*Math.Abs(det_F_k(k));
            for(int r = 0; r < M; r++) //Zeile
            {
                for(int s = 0; s < M; s++) //Spalte
                {
                   double[] x_1 = F_k_inverse_T(k, r);
                   double[] x_2 = F_k_inverse_T(k, s);
                    stiffness[r, s] = det*(x_1[0] * x_2[0] + x_1[1] * x_2[1]);
                }
            }
            return(stiffness);
        }

        public double[,] local_mass(int k)
        {
            double [,] mass=new double[M, M];
            double det = Math.Abs(det_F_k(k));
            for (int r = 0; r < M; r++) //Zeile
            {
                for (int s = 0; s < M; s++) //Spalte
                {
                    if (r == s)
                    {
                        mass[r, s] = det * (double)1 / (double)12;
                    }
                    else
                    {
                        mass[r,s]=det* (double)1 / (double)24;
                    }
                }
            }

            return (mass);
        }

      public double[,] global_mass()
        {
            double[,] C= new double[N,N];
            for (int k = 0; k < K; k++)
            {
                double[,] loc_mass=local_mass(k); 
                for(int i = 0; i < N; i++)
                {
                    for(int j = 0; j < N; j++)
                    {
                        for(int r = 0; r < M; r++)
                        {
                            for(int s = 0; s < M; s++)
                            {
                                if (A[r,k]==i & A[s, k] == j)
                                {
                                    C[i,j] += loc_mass[r, s];
                                 
                                }
                            }
                        }
                    }
                }  
            }
            return (C);
        }

       public double[,] global_stiffness()
        {
            double[,] C=new double[N,N];
            for (int k = 0; k < K; k++)
            {
                double[,] local_stiff = local_stiffness(k);
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int r = 0; r < M; r++)
                        {
                            for (int s = 0; s < M; s++)
                            {
                                if (A[r, k] == i & A[s, k] == j)
                                {
                                    C[i, j] += local_stiff[r, s];
                                }

                            }
                        }
                    }
                }
            }
            return (C);
        }

        public double[] stiffness_vec(double[] x)
        {
            double[]erg= new double[N];
            for(int i= 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    erg[i] += stiffness[i, j] * x[j];
                }
            }
            return (x);
        }


        public double[] local_load_vec(int k) //berechnet local load Vektor für Dreieck k (rechte Seite von Au=f)
        {
            double[] erg = new double[M];
            double det = Math.Abs(det_F_k(k));
            double x = (double)1 / (double)6;
            double[] x_1 = { x, x };
            double[] x_2 = { 4*x, x };
            double[] x_3 = { x, 4*x };
            erg[0] =det* x * (f(F_k(k, x_1)) * (1 - 2 * x) + f(F_k(k, x_2)) * (1 - 5 * x) + f(F_k(k, x_3)) * (1 - 5 * x)); //1-x-y
            erg[1] =  det*x * (f(F_k(k, x_1)) * x + f(F_k(k, x_2)) * 4*x + f(F_k(k, x_3)) * x); //x
            erg[2] = det* x * (f(F_k(k, x_1)) * x + f(F_k(k, x_2)) * x + f(F_k(k, x_3)) * 4*x);//y
            return (erg);
        }

        public double[] global_load_vec() //berechnet rechte Seite von Au=f
        {
            double[] erg = new double[N];
            for(int k= 0; k < K; k++)
            {
                double[] load_vec = local_load_vec(k);
                for(int i = 0; i< N; i++)
                {
                    for(int r = 0; r< M; r++)
                    {
                        if (A[r, k] == i)
                        {
                            erg[i] += load_vec[r];
                        }
                    }
                }
            }
            //Console.WriteLine("Counter ist: " + counter);
            return (erg);
        }

        public double[] F_k(int k, double[] x) //berechnet F_k(x), Transformation der Dreiecke
        {
            double[] erg=new double[2];
            erg[0] = (p(k, 1, 0) - p(k, 0, 0)) * x[0] + (p(k, 2, 0) - p(k, 0, 0)) * x[1] + p(k, 0, 0);
            erg[1] = (p(k, 1, 1) - p(k, 0, 1)) * x[0] + (p(k, 2, 1) - p(k, 0, 1)) * x[1] + p(k, 0, 1);
            return (erg);
        }

        public double det_F_k(int k) //berechnet det(F_k)
        {
            double det = (p(k, 1, 0) - p(k, 0, 0))*(p(k,2,1)-p(k,0,1))-(p(k,1,1)-p(k,0,1))*(p(k,2,0)-p(k,0,0));
            return (det);
        }

        public double[] F_k_inverse_T(int k, int r) //ist der Gradient von F_k^(-T), Matrix-Vektor Multiplikation, 
        {
            double[] erg = new double[2];
            double[,]B=new double [2,2];
            double det = det_F_k(k);
            B[0, 0] =1/det*( p(k, 2, 1) - p(k, 0, 1)); //Das ist F_k^(-1)
            B[0, 1] = 1/det*(p(k, 0, 0) - p(k, 2, 0));
            B[1, 0] = 1/det*(p(k, 0, 1) - p(k, 1, 1));
            B[1, 1] = 1/det*(p(k, 1, 0) - p(k, 0, 0));

            if (r == 0) //hier wird F_k^(-T)  grad\psi_r berechnet (-1,-1), (1,0), (0,1)
            {
                erg[0] = -B[0, 0] - B[1, 0];
                erg[1] = -B[0, 1] - B[1, 1];
            }
            if (r == 1)
            {
                erg[0] = B[0, 0];
                erg[1] = B[0, 1];
            }
            if (r == 2)
            {
                erg[0] = B[1, 0];
                erg[1] = B[1, 1];
            }
            return (erg);
        }

        public double[] P_V(double[] x) //Projector on V
        {
            double[] y = new double[N];
            for(int i = 0; i < N; i++)
            {
                y[i] = x[i];
            }
            for (int i = 0; i < D_index.Length; i++)
            {
                y[D_index[i]] = 0;
            }
            return(y);
        }

        public double [] P_D(double[]  y) //Projector on D
        {
            double[] erg = new double[N];
            double[] PV = P_V(y);
            for(int i=0; i < N; i++)
            {
                erg[i] =y[i]- PV[i];
            }
            
            return(erg);
   
        }


        public double[] global_stiffness_vec_boundary(double[] x)
        {
            double[] y_1 = P_V(stiffness_vec(P_V(x))); //P_V*A*P_V*x
            double[] y_2 = P_D(x); //P_D *x
            y_1=vec_add(y_1,y_2);
            return(y_1);
        }

        public double[] global_load_vec_boundary() //rechte Seite von Au=f
        {
            double[] u_D= new double[N];
            for(int i=0; i < N; i++)
            {
                u_D[i] = 0;
            }
            for(int i=0; i< D_index.Length; i++)
            {
                double[] y = { points[0, D_index[i]], points[1, D_index[i]] };
                u_D[D_index[i]] = f(y); 
            }
            double[] y_1 = P_V(vec_subtract(global_load_vec(), stiffness_vec(u_D)));
            y_1 = vec_add(y_1, u_D);
            return (y_1);
        }

        public double L2error(double[] x)
        {
            double erg = 0;
            double[,] A = global_mass();
            for(int i=0; i<N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    erg += x[j]* A[i, j] * x[j];
                }
            }
            return (erg);
        }


        public double[] vec_add(double[] x, double[] y)
        {
            for(int i = 0; i < x.Length; i++)
            {
                x[i] += y[i];  
            }
            return(x);
        }

        public double[] vec_subtract(double[] x, double[] y)
        {
            for (int i = 0; i < x.Length; i++)
            {
                x[i] -= y[i];
            }
            return (x);
        }




        /// <summary>
        /// Rechnet das Produkt von vektoren aus
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        private double produkt(double[] x, double[] y) //deine CG-Methode
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

            vector = global_stiffness_vec_boundary(xk);
            for (int j = 0; j < dimension; j++)
            {
                trk[j] = f[j] - vector[j];
                tdk[j] = trk[j];
            }
            do
            {
                vector = global_stiffness_vec_boundary(tdk);
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
                trk1 = global_stiffness_vec_boundary(xk);
                for (int j = 0; j < dimension; j++)
                {
                    vector[j] = f[j] - trk1[j];
                }
                fehler = produkt(vector, vector);

                i++;
            } while (i < Iteration && fehler > eps);

            return xk;

        }



















        }
    }

