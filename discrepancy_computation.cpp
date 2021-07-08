#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <stdlib.h>
#include <math.h>
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

// Computes the bound for discrepancy for specific m
inline double discrepancy(int p, int c, int d, int m)
{
   vector<double> val_set(p, 0);
   for (int i = 0; i < p; i++)
      val_set[i] = ((i + c * Pow(i, d, p)) % p)/(double)p;
   
   double msummation = 0;

   for (int h = 1; h <= m; h++)
   {  
      complex<double> Nsummation(0, 0);
      for (int i = 0; i < p; i++)
      {
         Nsummation = Nsummation + std::exp(2*pi*h*val_set[i]*ii);
      }
      double Nnorm = std::abs(Nsummation)/(double)p;
      msummation += (1/(double)h - 1/(double)(m + 1)) * Nnorm;
   }
   double disc = 6/(double)(m + 1) + 4/pi * msummation;
   return disc;
}

// Computes the minimum bound for discrepancy assuming it is a valley
double min_discrepancy(int p, int c, int d)
{
   int left = 6, right = 2*p;
   double discl, discr;
   while (left < right)
   {
      // if it causes integer overflow then it will run too long
      discl = discrepancy(p, c, d, (left + right)/2); 
      discr = discrepancy(p, c, d, (left + right)/2+1);
      if (discl < discr)
      {
         right = (left + right)/2;
      }
      else
         left = (left + right)/2 + 1;
   }
   return (discl < discr ? discl : discr);
}

void discrepancy_pc(int p)
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
            double min_disc = min_discrepancy(p, c, d);
            rec << "degree = " << d << "; c = " << c << "; min disc = " << min_disc << ";" << endl;
         }
      }
   }
   rec.close();
}

int main(int argc, char **argv)
{   
   if ( argc == 5 )
   {
      int a[ 4 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) m\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      cout << discrepancy(a[0], a[1], a[2], a[3]) << endl;
   }
   else if ( argc == 2 )
   {
      int a[ 1 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for list of discrepancy for cx^d+x over F_p:  <executable> p(prime)\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      discrepancy_pc(a[0]);
   }
   else
   {
      fprintf ( stderr, "Usage for discrepancy for cx^d+x over F_p:  <executable> p(prime) c(coefficient for x^d) d(degree) m\n" );  
      return -1;  
   }
}