#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#include "common.h"

using namespace std;


void functional_graph(int p, int a, int c, int d)
{
   vector<int> adj_mat(p, 0);
   fstream func_gr;
   func_gr.open(to_string(p) + "_" + to_string(a) + "_" + to_string(c) + "_" + to_string(d) + "_func_gr.csv", fstream::out);
   for (int i = 0; i < p; i++)
      adj_mat[i] = (a + i + c * Pow(i, d, p)) % p;
   
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
   else
   {
      fprintf ( stderr, "Usage for functional graph for cx^d+x+a over F_p:  <executable> p(prime) a(constant) c(coefficient for x^d) d(degree)\n" );  
      return -1;  
   }
}