
using System;
using System.Diagnostics.Metrics;
using System.Globalization;
using System.Runtime.CompilerServices;
using Test;
using System.IO;

int dim = 15;
int N = dim * dim; //Anzahl nodes in Triangulation
int K = 2 * (dim - 1) * (dim - 1); //Anzahl Dreiecke in Triangulation
int [,] A= new int [3,K]; //Matrix A speichert globale Index der Dreiecke
double[,] points= new double[2,N]; //speichert Koordinaten der Punkte der globalen Indices
int[] D_index = new int[4*dim-4]; //speichert die Dirichlet Indices;
for(int i=0; i<dim; i++)
{
    D_index[2*i] = i;
    D_index[2*i + 1] = dim * (dim - 1) + i;
}
for(int i=0; i<(dim-2); i++)
{
    D_index[2 * dim + 2*i] = dim + dim  * i;
    D_index[2 * dim + 2*i + 1] = 2 * dim - 1 + dim * i;
}

double f(double[] x)
{
    return  2 * ((1 - x[1]) * x[1] + (1 - x[0]) * x[0]);
}
double g(double[] x)
{
    return 0;
}

for(int i = 0; i < dim; i++) //Zeilenwert
{
    for(int j=0; j < dim; j++) //Spaltenwert
    {
        points[0,j * dim + i] = (double)i/(double)(dim - 1);
        points[1,j * dim + i] = (double)j / (double)(dim - 1);
    }
}

int counter = 0;
for(int i = 0; i< dim-1; i++) 
{
    for(int j = 0; j < 2 * (dim - 1);  j++)
    {
        A[0, counter] = dim * i + j % (dim-1);
        A[1, counter] = dim * i + j + 1;
        A[2, counter] = dim*i+ (dim+1) + j % (dim-1);
        counter += 1;
        //Console.WriteLine("counter"+counter);
        //Console.WriteLine("Index ist:" + A[0,counter-1]); 
    }
}



FEM_3012 T = new FEM_3012(A, points, D_index,f,g);



/*double[,] D = T.global_stiffness();
double[] a = new double[N];
for(int i=0; i<N; i++)
{
    a[i] = 3;
}
double [] C = T.F_k(0,a);
double[] d = T.F_k_inverse_T(0, 0);
double[] z = T.global_stiffness_vec_boundary(a);*/

/*for (int i = 0; i < D.GetLength(0); i++)
{
    for (int j = 0; j <D.GetLength(1); j++)
    {
        Console.Write("{0} ", D[i, j]);
    }
    Console.WriteLine();
}

for(int i = 0; i <z.Length; i++)
{
    Console.WriteLine("Eintrag :"+i+" ist " + z[i]);
}*/

double[] u = new double[N];
double[] zero = new double[N];
for (int i = 0; i < zero.Length; i++)
{
    zero[i] = 1;
    u[i] = 0;
}
double zer0 = 0;
Console.WriteLine("Berechnet CG jetzt");
u = T.CG_method(ref zer0, zero, T.P_V(T.global_load_vec()), ref dim);
Console.WriteLine("Vektor u von Au=f: ");
for (int i = 0; i < N; i++)
{
    Console.WriteLine(u[i]);
}

Console.WriteLine("det: "+ T.det_F_k(0));

StreamWriter sr = new StreamWriter("C:\\Users\\alexb\\Desktop\\FEM_project\\FEM_3012\\test.txt", false);
CultureInfo myCI = new CultureInfo("en-US", false);
for (int i=0; i < N; i++)
{
    sr.WriteLine(u[i].ToString()+ ";"+ points[0,i].ToString()+";"+ points[1,i].ToString(),myCI);
}
sr.Close();


