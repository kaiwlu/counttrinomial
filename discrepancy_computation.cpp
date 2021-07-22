#include <complex>
#include <iostream>
#include <fstream>
#include <set>
#include <stdlib.h>
#include <string>

#include "common.h"

using std::complex;
using std::fstream;
using std::vector;
using std::cout;
using std::endl;
using std::set;

const double pi = std::acos(-1);
const complex<double> ii(0, 1);

// Computes the bound for discrepancy of a sequence for specific m
inline double discrepancy(vector<double> &sequence, int m)
{
   
   double msummation = 0;

   for (int h = 1; h <= m; h++)
   {  
      complex<double> Nsummation(0, 0);
      for (size_t i = 0; i < sequence.size(); i++)
      {
         Nsummation = Nsummation + std::exp(2*pi*h*sequence[i]*ii);
      }
      double Nnorm = std::abs(Nsummation)/(double)(sequence.size());
      msummation += (1/(double)h - 1/(double)(m + 1)) * Nnorm;
   }
   double disc = 6/(double)(m + 1) + 4/pi * msummation;
   return disc;
}

// Computes the bound for discrepancy for cx^d+x
inline double binom_discrepancy(int p, int c, int d, int m, bool include_0)
{
   vector<double> val_set;
   if(include_0)
   {
      val_set.resize(p);
      for (int i = 0; i < p; i++)
         val_set[i] = ((i + c * Pow(i, d, p)) % p)/(double)p;
   }
   else
   {
      for (int i = 0; i < p; i++)
      {
         int eval = ((i + c * Pow(i, d, p)) % p);
         if (eval)
         {
            val_set.push_back((double)eval/(double)p);
         }
      }
   }

   return discrepancy(val_set, m);
}

// Computes the minimum bound for discrepancy assuming it is a valley
double min_discrepancy(int p, int c, int d, bool include_0)
{
   int left = 6, right = 2*p;
   double discl, discr;
   while (left < right)
   {
      // if it causes integer overflow then it will run too long
      discl = binom_discrepancy(p, c, d, (left + right)/2, include_0); 
      discr = binom_discrepancy(p, c, d, (left + right)/2+1, include_0);
      if (discl < discr)
      {
         right = (left + right)/2;
      }
      else
         left = (left + right)/2 + 1;
   }
   return (discl < discr ? discl : discr);
}

void discrepancy_pc(int p, bool include_0)
{
   fstream rec;
   rec.open("min_disc_" + std::to_string(p) + ".txt", fstream::out);
   set<int> coeffc;
   unsigned int coeffsize = (int)pow(p, 0.25) + 2;
   while (coeffc.size() < coeffsize)
   {
      coeffc.insert(rand() % (p - 2) + 1); // c ranging from 1 to P - 2
   }
   set<int> expd;
   while (expd.size() < (unsigned int)p/5)
   {
      expd.insert(rand() % (p - 4) + 3); // d ranging from 3 to P - 2
   }
   for (int d : expd)
   {
      if (Coprime(d, p - 1))
      {
         for (int c : coeffc)
         {
            double min_disc = min_discrepancy(p, c, d, include_0);
            rec << "degree = " << d << "; c = " << c << "; min disc = " << min_disc << ";" << endl;
         }
      }
   }
   rec.close();
}

void discrepancy_graph(int N, bool include_0)
{
   vector<int> prime_list = Primes(0, N);
   fstream disc_gr;
   if (include_0)
   {
      disc_gr.open("discrepancy_data.csv",fstream::out);
   }
   else
   {
      disc_gr.open("discrepancy_data_n0.csv",fstream::out);
   }
   
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      for (int d = (int)pow(prime_list[i], 0.6); d < prime_list[i] - 1; d++)
      {
         if (!Coprime(d, prime_list[i]-1))
            continue;

         for (int c = 1; c < prime_list[i]; c++)
         {
            disc_gr << prime_list[i] << "," << c << "," << d << "," << min_discrepancy(prime_list[i], c, d, include_0) << "\n";
         }
      }
      disc_gr.flush();
   }  
   disc_gr.close();
}

int main(int argc, char **argv)
{   
   if ( argc == 6 )
   {
      int a[ 4 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) m <include 0 ? 'Y' : 'N'>\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      cout << binom_discrepancy(a[0], a[1], a[2], a[3], (*argv[5] == 'Y')) << endl;
   }
   if ( argc == 5 )
   {
      int a[ 3 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for min discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) <include 0 ? 'Y' : 'N'>\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      cout << min_discrepancy(a[0], a[1], a[2], (*argv[4] == 'Y')) << endl;
   }
   else if ( argc == 4 )
   {
      int N;  

      char* p_end;  
      const long element = std::strtol( argv[ 2 ], &p_end, 10 );  
      if ( argv[ 2 ] == p_end || argv[ 1 ][ 0 ] != 'G')
      {
         fprintf ( stderr, "Usage for discrepancy data for cx^d+x over F_p for p up to N:  <executable> 'G' N <include 0 ? 'Y' : 'N'>\n" );  
         return -1;  
      }
      N = element;  
      discrepancy_graph(N, (*argv[3] == 'Y'));
   }
   else if ( argc == 3 )
   {
      int a[ 1 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for list of discrepancy for cx^d+x over F_p:  <executable> p(prime) <include 0 ? 'Y' : 'N'>\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      discrepancy_pc(a[0], (*argv[3] == 'Y'));
   }
   else
   {
      fprintf ( stderr, "Usage for discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) m <include 0 ? 'Y' : 'N'>\nUsage for min discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) <include 0 ? 'Y' : 'N'>\nUsage for discrepancy data for cx^d+x over F_p for p up to N:  <executable> 'G' N <include 0 ? 'Y' : 'N'>\nUsage for list of discrepancy for cx^d+x over F_p:  <executable> p(prime) <include 0 ? 'Y' : 'N'>\n" );  
      return -1;  
   }
}