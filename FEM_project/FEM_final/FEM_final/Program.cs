using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;


namespace FEM_final
{

 /*
  * Es feht noch deine CG-Methode
  * Du musst das Gleichungssystem Au=f lösen, wobei
  * f durch D.global_load_vec() gegeben ist (der Vektor wird auch unten ausgegeben)
  * Die Matrixvektormultiplikation Au gegeben durch stiffness_vec(x)
  * 
  * 
  * 
  */
    internal class Program
    {
        static void Main(string[] args)
        {

            int N = 5; //Anzahl Knoten in einer Reihe; hier wird eine symmetrische Triangulation erstellt wie auf dem Blatt eingezeichnet
            double[] x_1= new double [12 * (N - 1) * (N - 1)];
            int k = 0;
            for(int i=0; i < N - 1; i++) //hier werden alle Dreiecke in den Vektor x_1 geschrieben; ein Dreieck wird durch 6 Punkte definiert
            {
                for(int j=0; j < N - 1; j++)
                {
                    x_1[12 * k] = i;
                    x_1[12 * k + 1] = j;
                    x_1[12 * k + 2] = i + 1;
                    x_1[12 * k + 3] = j;
                    x_1[12 * k + 4]= i + 1;
                    x_1[12 * k + 5]= j + 1;
                    x_1[12 * k + 6] = i;
                    x_1[12 * k + 7] = j;
                    x_1[12 * k + 8]= i + 1;
                    x_1[12 * k + 9]=j + 1;
                    x_1[12 * k + 10] = i;
                    x_1[12 * k + 11] = j + 1;
                    k++;
                }
            }
            double[] DB= new double[8*(N-2)]; // hier werden die Punkte mit Dirichlet-boundary festgelegt; in unserem Fall alle Randpunkte außer die 4 Eckpunkte des Quadrats
            int l = 0;
            for(int i = 1; i < N - 1; i++)
            {
                DB[8*l] = i;
                DB[8*l + 1] = 0;
                DB[8 * l + 2] = 0;
                DB[8 * l + 3] = i;
                DB[8 * l + 4] = i;
                DB[8 * l + 5] = N-1;
                DB[8 * l + 6] = N-1;
                DB[8 * l + 7] = i;
                l++;
            }
            double[] NB = { 0, 0, 0, N-1, N-1, 0, N-1, N-1 }; //Neunmann-boundary sind die 4 Eckpunkte des Quadrats
            FEM D = new FEM(x_1, 2,NB,DB);
            //D.print_points();
            D.print_global_indeces();
            D.print_set_i_N();
            D.print_set_D();
            D.number_triangles();

            Console.WriteLine("Vektor f von Au=f: ");
            for(int i = 0; i < N*N; i++)
            {
                Console.WriteLine(D.global_load_vec()[i]);
            }

 
        }
    
    }
}

