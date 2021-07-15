#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#include "common.h"

using namespace std;
static int largest_p;

// PLEASE INITIALIZE largest_p before calling this function for the first time
int primitive_root(int p)
{
   static vector<int> prime_list = Primes(1, largest_p);
   vector<int> powers;
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      if (prime_list[i] <= p)
      {
         if ((p - 1) % prime_list[i] == 0)
         {
            powers.push_back((p - 1)/prime_list[i]);
         }
      }
      else
         break;
   }

   int gent = -1;

   for (int i = 2; (i < p) && (gent < 0); i++)
   {
      gent = i;
      for (size_t j = 0; j < powers.size(); j++)
      {
         if (Pow(i, powers[j], p) == 1)
         {
            gent = -1;
            break;
         }
      }
   }
   
   return gent;
}


void functional_graph(int p, int a, int c, int d)
{
   vector<int> adj_mat(p, 0);
   fstream func_gr;
   func_gr.open(to_string(p) + "_" + to_string(a) + "_" + to_string(c) + "_" + to_string(d) + "_func_gr.csv", fstream::out);
   for (int i = 0; i < p; i++)
      adj_mat[i] = (a + i + Mult(c, Pow(i, d, p), p)) % p;
   
   for (int i = 0; i < p; i++)
   {
      for (int j = 0; j < p; j++)
      {
         if(adj_mat[i] == j)
            func_gr << "1" << (j == p - 1 ? "\n" : ", ");
         else 
            func_gr << "0" << (j == p - 1 ? "\n" : ", ");
      }
   }
   func_gr.close();
}

void exp_functional_graph(int p)
{
   vector<int> adj_mat(p, 0);
   fstream func_gr;
   func_gr.open("exp_" + to_string(p) + "_func_gr.csv", fstream::out);
   largest_p = p;
   int g = primitive_root(p);
   if (g == -1)
   {
      func_gr.close();
      fprintf ( stderr, "primitive root is wrong\n" );
      return;
   }
   
   for (int i = 0; i < p; i++)
      adj_mat[i] = Pow(g, i, p) % p;
   
   for (int i = 0; i < p; i++)
   {
      for (int j = 0; j < p; j++)
      {
         if(adj_mat[i] == j)
            func_gr << "1" << (j == p - 1 ? "\n" : ", ");
         else 
            func_gr << "0" << (j == p - 1 ? "\n" : ", ");
      }
   }
   func_gr.close();
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
            fprintf ( stderr, "Usage for functional graph for cx^d+x+a over F_p:  <executable> p(prime) a(constant) c(coefficient for x^d) d(degree)\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      functional_graph(a[0], a[1], a[2], a[3]);
   }
   else if ( argc == 2 )
   {
      int p;
      char* p_end;  
      const long element = std::strtol( argv[ 1 ], &p_end, 10 );  
      if ( argv[ 1 ] == p_end )
      {
         fprintf ( stderr, "Usage for functional graph for g^x over F_p:  <executable> p(prime)\n" );  
         return -1;  
      }
      p = element;
      exp_functional_graph(p);
   }
   else
   {
      fprintf ( stderr, "Usage for functional graph for cx^d+x+a over F_p:  <executable> p(prime) a(constant) c(coefficient for x^d) d(degree)\n" );  
      return -1;  
   }
}