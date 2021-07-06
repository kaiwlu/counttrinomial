#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

// Integer powers (x^exp) mod p for p up to 46340, above that it might overflow
int Pow(int x, int exp, int p)
{
  if (exp == 0) return 1;
  if (exp == 1) return x;
  
  int tmp = Pow(x, exp/2, p);
  if (exp % 2 == 0) return (tmp * tmp) % p;
  else return (((x * tmp) % p) * tmp) % p;
}

class Mat
{
private:
   int *elements;
   int width, height;

public:
   Mat(int w, int h, int init)
   {
      elements = new int[ width * height ];
      width = w;
      height = h;
      for (int i = 0; i < width * height; i++)
      {
         elements[i] = init;
      }
   }
   ~Mat()
   {
      delete[] elements;
   }
   int& elt(int x, int y)
   {
      return elements[ x + width * y ];
   }
};

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