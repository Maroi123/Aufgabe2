using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Globalization;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace WissRech
{
    public class CG
    {
        public double get_fehler()
        {
            return fehler;
        }
        public double[] get_solution()
        {
            return solution;
        }
        public void set_b_x_(double[] b, double[] x)
        {
            if (b.Length != x.Length)
            {
                throw new InvalidOperationException("Vektoren haben nicht dieselbe Dimension");
            }
            else
            {
                this.b = b;
                this.xk = x;
                this.dimension = b.Length;
            }
        }
        public void set_verfahren(string Methode)
        {
            this.Methode = Methode;
        }
        public void set_Anzahl(int N)
        {
            this.N = N;
        }
        public void run()
        {
            switch (Methode)  //wähle die Methode aus
            {
                default:   //Erklärung: Switch methode funktioniert fast wie ein if nur mit dem unterschied das in den richtigen fall gesprungen wird bsp: du hast 100 fälle du könntest immer if(i==1){code} elseif(i==2){code} der computer prüft dann jedes mal ist i=1? nein dann ist i=2? nein usw.. switch ist einfach effektiver
                    break;
                case "Hilbert":
                    CG_method(hilbert);  //führe das CG verfahren mit Hilbert aus
                    break;
                case "Identity":
                    CG_method(identity); //führe das CG verfahren mit der id aus
                    break;
                case "FE":
                    CG_method(FE); //führe das CG verfahren mit der FE aus
                    break;



            }
        }
        public void Ausgabe()
        {
            Console.WriteLine("Der Lösungsvektor ist: ");
            for (int j = 0; j < dimension; j++) //gebe den lösungsvektor aus
            {
                solution[j] = xk[j];
                Console.WriteLine(solution[j]);
            }
            Console.WriteLine("Der fehler beträgt ca: " + fehler); //und den fehler
        }
        public void set_fehler(double eps)
        {
            this.eps = eps;
        }
       
  
        double[] b, solution, xk1, xk, tdk, tdk1, trk, trk1; //b = rechte seite, solution ist der Lösungsvektor xk1= xk+1, xk=xk,zk =z vektor ,tdk=dk,tdk1=tdk+1,trk=rk, trk1=trk+1;
        double eps, fehler, ak, bk; //eps ist die genauigkeit, fehler ist der aktuelle fehler, ak ist alpha k, bk ist beta k
        int N, dimension;    //N ist die anzahl der Iterationen, dimension ist die Dimension der Vektoren
        string Methode; //Welche Methode verwendet werden soll
        /// <summary>
        /// Hier wird ein Objekt von CG konstruiert
        /// </summary>
        /// <param name="b">Vektor der rechten seite Ax=b</param>
        /// <param name="xk">startwert</param>
        /// <param name="eps">Fehlergenauigkeit</param>
        /// <param name="N">Iterationen</param>
        /// <param name="Methode">Welche Methode verwendet werden soll</param>
        public CG(double[] b, double[] xk, double eps, int N, string Methode)
        {
            if (b.Length != xk.Length)
            {
                throw new InvalidOperationException("dimension not fitting");  //werfe einen fehler wenn die dimensionen nicht passen. 
            }
            else
            {
                this.dimension = b.Length; //setzte die dimension
                objektifizierung(); // erstelle die arrays
                this.N = N; //setze die anzahl der durchläufe
                this.b = b; //setze die rechte seite
                this.xk = xk; //setze den startvektor
                this.eps = eps;  //setze den maximalen fehler den man haben will, dies ist ein oder, entweder der fehler wird erreicht oder die anzahl der durchläufe
                this.Methode = Methode; //setze die methode
                

                
                for (int i = 0; i < dimension; i++) //initialisiere die array's mit 0 dieser schritt ist semi unnötig
                {
                    this.solution[i] = 0;
                    this.xk1[i] = 0;
                    
                    this.tdk[i] = 0;
                    this.tdk1[i] = 0;
                    this.trk[i] = 0;
                    this.trk1[i] = 0;
                    
                 


                }

                
               
            }

            

        }
        
        /// <summary>
        /// objekte von den Vektoren erstellen. 
        /// </summary>
        private void objektifizierung()
        {
            this.b =new double[dimension];
            this.solution = new double[dimension];
            this.xk1 = new double[dimension];
            this.xk = new double[dimension];
       
            this.tdk = new double[dimension];
            this.tdk1 = new double[dimension];
            this.trk = new double[dimension];
            this.trk1 = new double[dimension];
            
         


        }
        /// <summary>
        /// implementiert die CG methode
        /// </summary>
        /// <param name="matrix_vec">Welche Matrix verwendet werden soll</param>
        private void CG_method(Func<double[], double[]> matrix_vec)
        {

            int i = 0;
            double[] vector = new double[dimension];
            vector = matrix_vec(xk);
            for(int j = 0; j < dimension; j++)
            {
                trk[j] = b[j] - vector[j];
                tdk[j] = trk[j];
            }
            do
            {
                vector = matrix_vec(tdk);
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
                trk1 = matrix_vec(xk);
                for (int j = 0; j < dimension; j++)
                {
                    vector[j] = b[j] - trk1[j];
                }
                fehler = produkt(vector, vector);

                i++;
            } while (i < N && fehler > eps);
           // Console.WriteLine("Das ergebnis:Ax=");
            //vector = matrix_vec(xk);
            //for(int j=0;j<dimension;j++)
            //{
              //  Console.WriteLine(vector[j]);
            //}
            
        }
        /// <summary>
        /// berechnet transpose(x)*y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private double produkt(double[] x, double[] y) 
        {
            double outp = 0;
            if(x.Length != y.Length)
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
        /// Multiplikation mit id
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double[] identity(double[] x)
        {
            return x;
        }
        /// <summary>
        /// Multiplikation mit der Hilbertmatrix
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double[] hilbert(double[] x)
        {
            double sum = 0;
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                sum = 0;
                for (int j = 0; j < x.Length; j++)
                {
                    sum += x[j] * 1 / (i + j + 1);
                }
                y[i] = sum;
            }
            return y;
        }


        /// <summary>
        /// Multiplikation mit der Finiten Elementen Matrix
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double[] FE(double[] x)
        {
            double[] y = new double[x.Length];
            for (int i = 1; i < x.Length - 1; i++)
            {
                y[i] = -x[i - 1] + 2 * x[i] - x[i + 1];
            }
            y[0] = 2 * x[0] - x[1];
            y[x.Length - 1] = -x[x.Length - 2] + 2 * x[x.Length - 1];
            return y;
        }
    }
}
