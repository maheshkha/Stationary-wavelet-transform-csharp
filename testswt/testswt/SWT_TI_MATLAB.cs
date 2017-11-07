//#define DEBUG
#define MYTEST
using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;


namespace testswt
{
    class Program
    {

        static double[] lifting_step(double[] x, double[] h, int dir)
        {
            if (dir == 1)
            {
                int N = x.Length;
                double[] evenList = new double[N / 2];
                double[] oddList = new double[N / 2];
                double[] tmpevenList = new double[N / 2];
                double[] tmpoddList = new double[N / 2];
                int temp = 0, temp1 = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (i % 2 == 0)
                    {
                        evenList[temp] = x[i];
                        temp++;
                    }
                    else
                    {
                        oddList[temp1] = x[i];
                        temp1++;
                    }
                }

                for (int m1 = 0; m1 < 2; m1++)
                {
                    for (int i = 0; i < tmpevenList.Length; i++)
                    {
                        if (i < tmpevenList.Length - 1)
                        {
                            tmpevenList[i] = evenList[i + 1];
                        }
                        else
                        {
                            tmpevenList[tmpevenList.Length - 1] = evenList[tmpevenList.Length - 1];
                        }
                    }

                    for (int i = 0; i < tmpevenList.Length; i++)
                    {
                        oddList[i] = oddList[i] - h[2 * m1] * (evenList[i] + tmpevenList[i]);
                    }

                    for (int i2 = 0; i2 < oddList.Length; i2++)
                    {
                        if (i2 > 0)
                        {
                            tmpoddList[i2] = oddList[i2 - 1];
                        }
                        else
                        {
                            tmpoddList[i2] = oddList[0];
                        }
                    }

                    for (int i = 0; i < evenList.Length; i++)
                    {
                        evenList[i] = evenList[i] + h[2 * m1 + 1] * (oddList[i] + tmpoddList[i]);
                    }
                }

                for (int i = 0; i < evenList.Length; i++)
                {
                    evenList[i] = evenList[i] * h[h.Length - 1];
                    oddList[i] = oddList[i] / h[h.Length - 1];
                }

                for (int i = 0; i < x.Length / 2; i++)
                {
                    x[i] = evenList[i];
                    x[i + (x.Length / 2)] = oddList[i];
                }

                return x;
            }
            else
            {
                int N = x.Length;
                double[] evenList = new double[N / 2];
                double[] oddList = new double[N / 2];
                double[] tmpevenList = new double[N / 2];
                double[] tmpoddList = new double[N / 2];
                int temp = 0;
                for (int i = 0; i < x.Length / 2; i++)
                {
                    evenList[temp] = x[i + (x.Length / 2)] * h[h.Length - 1];
                    oddList[temp] = x[i] / h[h.Length - 1];
                    temp++;
                }

                for (int m1 = 2; m1 >= 1; m1--)
                {
                    for (int i = 0; i < tmpevenList.Length; i++)
                    {
                        if (i > 0)
                        {
                            tmpevenList[i] = evenList[i - 1];
                        }
                        else
                        {
                            tmpevenList[i] = evenList[0];
                        }
                    }
                    for (int i = 0; i < oddList.Length; i++)
                    {
                        oddList[i] = oddList[i] - h[2 * m1 - 1] * (evenList[i] + tmpevenList[i]);
                    }

                    for (int i2 = 0; i2 < oddList.Length; i2++)
                    {
                        if (i2 < oddList.Length - 1)
                        {
                            tmpoddList[i2] = oddList[i2 + 1];
                        }
                        else
                        {
                            tmpoddList[tmpoddList.Length - 1] = oddList[tmpoddList.Length - 1];
                        }
                    }
                    for (int i = 0; i < evenList.Length; i++)
                    {
                        evenList[i] = evenList[i] + h[2 * m1 - 2] * (oddList[i] + tmpoddList[i]);
                    }
                }

                for (int i = 0; i < x.Length; i++)
                {
                    if (i % 2 == 0)
                    {
                        x[i] = oddList[i / 2];
                    }
                    else
                    {
                        x[i] = evenList[i / 2];
                    }
                }

                return x;
            }
        }

        static double[, ,] lifting_step_ti(double[, ,] x, double[] h, int dir, int dist)
        {
            int m = (h.Length - 1) / 2;

            int N = x.GetLength(0);
            int[] s1 = new int[N];
            int[] s2 = new int[N];
            double[, ,] outputList = new double[N, 1, 2];
            double[] evenList = new double[N];
            double[] oddList = new double[N];

            int temp = 0, temp1 = 0;

            for (int i = 0; i < N; i++)
            {
                s1[i] = i + 1 + dist;
                s2[i] = i + 1 - dist;
                // boundary conditions
                if (s1[i] > N)
                {
                    s1[i] = 2 * N - s1[i];
                }
                if (s2[i] > N)
                {
                    s2[i] = 2 * N - s2[i];
                }
                if (s1[i] < 1)
                {
                    s1[i] = 2 - s1[i];
                }
                if (s2[i] < 1)
                {
                    s2[i] = 2 - s2[i];
                }
                s1[i] = s1[i] - 1;
                s2[i] = s2[i] - 1;
            }

            if (dir == 1)
            {
                for (int i = 0; i < N; i++)
                {
                    evenList[temp] = x[i, 0, 0];
                    oddList[temp1] = x[i, 0, 0];
                    temp++;
                    temp1++;
                }

                for (int m1 = 0; m1 < 2; m1++)
                {
                    for (int i = 0; i < evenList.Length; i++)
                    {
                        oddList[i] = oddList[i] - h[2 * m1] * (evenList[s1[i]] + evenList[s2[i]]);
                    }
                    for (int i = 0; i < evenList.Length; i++)
                    {
                        evenList[i] = evenList[i] + h[2 * m1 + 1] * (oddList[s1[i]] + oddList[s2[i]]);
                    }

                }

                for (int i = 0; i < evenList.Length; i++)
                {
                    evenList[i] = evenList[i] * h[h.Length - 1];
                    oddList[i] = oddList[i] / h[h.Length - 1];
                }

                for (int i = 0; i < evenList.Length; i++)
                {
                    outputList[i, 0, 0] = evenList[i];
                    outputList[i, 0, 1] = oddList[i];
                }
            }
            else
            {
                for (int i = 0; i < N; i++)
                {
                    evenList[temp] = x[i, 0, 1] * h[h.Length - 1];
                    oddList[temp] = x[i, 0, 0] / h[h.Length - 1];
                    temp++;
                }

                for (int m1 = 2; m1 >= 1; m1--)
                {
                    for (int i = 0; i < evenList.Length; i++)
                    {
                        oddList[i] = oddList[i] - h[2 * m1 - 1] * (evenList[s1[i]] + evenList[s2[i]]);
                    }
                    for (int i = 0; i < evenList.Length; i++)
                    {
                        evenList[i] = evenList[i] + h[2 * m1 - 2] * (oddList[s1[i]] + oddList[s2[i]]);
                    }

                }

                for (int i = 0; i < evenList.Length; i++)
                {
                    outputList[i, 0, 0] = (evenList[i] + oddList[i]) / 2;
                }

            }
            return outputList;
        }

        static double[, ,] perform_wavelet_transf(double[, ,] x, int Jmin, int dir)
        {
            // filter coefficients
            double[] h = new double[5];
            h[0] = 1.586134342;
            h[1] = -.05298011854;
            h[2] = -.8829110762;
            h[3] = .4435068522;
            h[4] = 1.149604398;

            int levelmin = Jmin;

            // calculate levels
            int n1 = x.GetLength(0);
            int levelmax = (int)Math.Log(n1, 2) - 1;
            int nJ = levelmax - levelmin + 1;
            double[, ,] dataDouble_rep = new double[n1, 1, nJ + 1];
            double[, ,] dataDouble_out = new double[n1, 1, nJ + 1];

            if (dir == 1)
            {
                for (int j = 0; j <= nJ; j++)
                {
                    for (int i = 0; i < n1; i++)
                    {
                        dataDouble_rep[i, 0, j] = x[i, 0, 0];
                    }
                }
                for (int j = levelmax; j >= levelmin; j--)
                {
                    int levelsize = (int)Math.Pow(2, j + 1);
                    int dist = (int)Math.Pow(2, levelmax - j);
                    int ind = j - levelmin + 1;

                    // Calling the lifting forward transform code
                    dataDouble_out = Program.lifting_step_ti(dataDouble_rep, h, dir, dist);

                    // Copy intermedicate level outputs
                    for (int i1 = 0; i1 < n1; i1++)
                    {
                        dataDouble_rep[i1, 0, 0] = dataDouble_out[i1, 0, 0];
                        dataDouble_rep[i1, 0, ind] = dataDouble_out[i1, 0, 1];
                    }

                }
            }
            else
            {
                for (int j = 0; j <= nJ; j++)
                {
                    for (int i = 0; i < n1; i++)
                    {
                        dataDouble_rep[i, 0, j] = x[i, 0, j];
                    }
                }
                for (int j = levelmin; j <= levelmax; j++)
                {
                    int levelsize = (int)Math.Pow(2, j + 1);
                    int dist = (int)Math.Pow(2, levelmax - j);
                    int ind = j - levelmin + 1;
                    // Copy intermedicate level outputs
                    for (int i1 = 0; i1 < n1; i1++)
                    {
                        dataDouble_out[i1, 0, 0] = dataDouble_rep[i1, 0, 0];
                        dataDouble_out[i1, 0, 1] = dataDouble_rep[i1, 0, ind];
                    }

                    // Calling the lifting forward transform code
                    dataDouble_out = Program.lifting_step_ti(dataDouble_out, h, dir, dist);

                    // Copy intermedicate level outputs
                    for (int i1 = 0; i1 < n1; i1++)
                    {
                        dataDouble_rep[i1, 0, 0] = dataDouble_out[i1, 0, 0];
                    }
                }

            }
            return dataDouble_rep;
        }

        static void Main(string[] args)
        {

            int counter = 0;
            string line;
            int N = 512;
            double[] dataDouble = new double[N];
            double[] dataOut = new double[N];
            double[, ,] dataDouble_rep1 = new double[N, 1, 9];
            double[, ,] dataDouble_rep2 = new double[N, 1, 9];


            Stopwatch sw = new Stopwatch();

            // Read the file and display it line by line.
            System.IO.StreamReader rfile =
                new System.IO.StreamReader("test512.txt");
            System.IO.StreamWriter wfile =
                        new System.IO.StreamWriter("fwd_op_test512.txt");
            System.IO.StreamWriter wfile1 =
                         new System.IO.StreamWriter("inv_op_test512.txt");

            while ((line = rfile.ReadLine()) != null)
            {
                dataDouble[counter] = double.Parse(line);
                counter++;
            }
            for (int j = 0; j < 9; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    dataDouble_rep1[i, 0, j] = dataDouble[i];
                }
            }

            sw.Start();
            dataDouble_rep2 = Program.perform_wavelet_transf(dataDouble_rep1, 1, 1);
            dataDouble_rep1 = Program.perform_wavelet_transf(dataDouble_rep2, 1, -1);
            // Copy intermedicate level outputs
            for (int i1 = 0; i1 < N; i1++)
            {
                dataOut[i1] = dataDouble_rep1[i1, 0, 0];
            }

            foreach (double i in dataOut)
            {
                //Console.WriteLine(i.ToString());
                wfile1.WriteLine(i.ToString());
            }

            sw.Stop();
            Console.WriteLine("\nTime taken: {0}ms", sw.Elapsed.TotalMilliseconds);

            rfile.Close();
            wfile.Close();
            wfile1.Close();

            // Suspend the screen.
            Console.ReadLine();
        }
    }
}
