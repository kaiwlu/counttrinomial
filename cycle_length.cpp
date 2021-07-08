#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <set>
#include "common.h"

using namespace std;

struct Cycle_data
{
   int roots = 0; // number of roots
   int cycle_num = 0; // number of cycles
   int cycle_len = 0; // number of periodic points
   int cycle_min = 2147483647; // length of shortest cycle
   int cycle_max = 0; // length of longest cycle
};

// lists primes from min up to n using sieve of Eratosthenes
vector<int> Primes(int min, int n)
{
   vector<bool> numbers(n,true);
   vector<int> prime_list;
   int sqrtn = (int)sqrt(n);
   for (int i = 2; i < n; i++)
   {
      if (numbers[i])
      {
         if (i >= min)
            prime_list.push_back(i);
         if (i <= sqrtn)
            for (int j = 2; i*j < n; j++)
               numbers[i*j] = false;
      } 
   }
   return prime_list;
}

// Returns cycle data for 
// Trinomials of the form f(x)=a+x+cx^d over F_p
Cycle_data Trinomial_cycle(int p, int a, int c, int d)
{
   Cycle_data tri_cyc;

   // val is a p*3 2d array
   // val[i][j] is then rewritten as val[i*3+j]
   int *val = new int[p*3];
   // val[][0] being -1 means not visited
   // val[][1] is the position in cycle
   // val[][2] is i
   for (int i = 0; i < p*3; i++)
      val[i] = -1;
   
   for (int i = 0; i < p; i++)
   {
      int j = i;
      int pos = 0;
      while (val[j*3+0] == -1)
      {
         pos++;
         val[j*3+0] = (a + j + c * Pow(j, d, p)) % p;
         val[j*3+1] = pos;
         val[j*3+2] = i;
         j = val[j*3+0];
         if (j == 0)
            tri_cyc.roots++;
      }

      // if the cycle is detected in this iteration, record its length
      if (val[j*3+2] == i)
      {
         int cyc_len = pos - val[j*3+1] + 1;
         tri_cyc.cycle_len += cyc_len;
         tri_cyc.cycle_num++;
         if (cyc_len > tri_cyc.cycle_max)
            tri_cyc.cycle_max = cyc_len;
         if (cyc_len < tri_cyc.cycle_min)
            tri_cyc.cycle_min = cyc_len;
      }
   }
   delete [] val;
   return tri_cyc;
}

// Records the average number of elements in cycles(cyc_len), average number
// of cycles(cyc_num), average length of smallest cycles(cyc_min), and 
// average length of largest cycles(cyc_max) for trinomials in the form 
// f(x)=a+x+cx^d for every degree over F_p, Min <= p < N
// csv line format: prime, degree, data\n
void Count_cycle_average(int Min, int N)
{
   vector<int> prime_list = Primes(Min, N);
   fstream roots, num, len, min, max;
   num.open("avg_cyc_num.csv",fstream::out);
   len.open("avg_cyc_len.csv",fstream::out);
   min.open("avg_cyc_min.csv",fstream::out);
   max.open("avg_cyc_max.csv",fstream::out);
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      double count = prime_list[i] * (prime_list[i] - 1);
      for (int d = 2; d < prime_list[i] - 1; d++)
      {
         int length_sum = 0;
         int min_sum = 0;
         int max_sum = 0;
         int num_sum = 0;
         for (int a = 0; a < prime_list[i]; a++)
            for (int c = 1; c < prime_list[i]; c++)
            {
               Cycle_data tri_data = Trinomial_cycle(prime_list[i], a, c, d);
               length_sum += tri_data.cycle_len;
               min_sum += tri_data.cycle_min;
               max_sum += tri_data.cycle_max;
               num_sum += tri_data.cycle_num;
            }
         len << prime_list[i] << ',' << d << ',' << length_sum/count << '\n';
         min << prime_list[i] << ',' << d << ',' << min_sum/count << '\n';
         max << prime_list[i] << ',' << d << ',' << max_sum/count << '\n';
         num << prime_list[i] << ',' << d << ',' << num_sum/count << '\n';
      }
   }
   roots.close();
   num.close();
   len.close();
   min.close();
   max.close();
}

// Records the Cycle_data for trinomials in the form 
// f(x)=a+x+cx^d for every degree d coprime to p-1 over F_p, n <= p < N
// csv line format: prime, degree, a, c, num_roots, num_periodic_pts, min_cycle_len, max_cycle_len, num_cycles\n
void Count_coprime_cycle(int n, int N)
{
   vector<int> prime_list = Primes(n, N);
   fstream cycdata;
   cycdata.open("ind_cyc_data.csv",fstream::out);
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      for (int d = (int)pow(prime_list[i], 0.6); d < prime_list[i] - 1; d++)
      {
         if (!Coprime(d, prime_list[i]-1))
            continue;

         for (int a = 0; a < prime_list[i]; a++)
            for (int c = 1; c < prime_list[i]; c++)
            {
               Cycle_data tri_data = Trinomial_cycle(prime_list[i], a, c, d);
               cycdata << prime_list[i] << ',' << d << ',' << a << ',' << c << ',' << tri_data.roots << ',' << tri_data.cycle_len << ',' << tri_data.cycle_min << ',' << tri_data.cycle_max << ',' << tri_data.cycle_num << endl;
            }

      }
      cycdata.flush();
   }  
   cycdata.close();
}

// Input p the prime you want cycle info over
// Records the individual cycle lengths for pow(p, 0.7) randomly selected trinomials of the form 
// f(x)=a+x+cx^d for every degree d >= pow(p, 0.6) coprime to p-1 over F_p
// Assuming a = 1 or 0, every cycle appears exactly once
// csv line format: d, a, c, cycle_len, cycle_len, ..., cycle_len\n
void Count_large_prime_cycle(int P)
{
   fstream cycdata;
   cycdata.open(to_string(P) + "_cyc_data.csv",fstream::out);

   set<int> coeffc;
   unsigned int coeffsize = (int)pow(P, 0.7) + 2;
   while (coeffc.size() < coeffsize)
   {
      coeffc.insert(rand() % (P - 2) + 1); // c ranging from 1 to P - 2
   }
   
   // val is a p*3 2d array
   // val[i][j] is then rewritten as val[i*3+j]
   vector<int> val(P*3, -1);

   for (int d = (int)pow(P, 0.6); d < P - 1; d++)
   {
      if (!Coprime(d, P-1))
         continue;
      for (int a : {0, 1})
         for (int c : coeffc)
         {
            cycdata << d << ',' << a << ',' << c;

            // val[][0] being -1 means not visited
            // val[][1] is the position in iteration, needed to calculate cycle length
            // val[][2] is i, needed to determine if cycle was detected in this iteration
            for (int i = 0; i < P*3; i++)
               val[i] = -1;
            
            int j;

            for (int i = 0; i < P; i++)
            {
               if (val[i*3+0] != -1)
                  continue;
               
               j = i;
               int pos = 0;
               while (val[j*3+0] == -1)
               {
                  pos++;
                  val[j*3+0] = (a + j + c * Pow(j, d, P)) % P;
                  val[j*3+1] = pos;
                  val[j*3+2] = i;
                  j = val[j*3+0];
               }
               // if the cycle is detected in this iteration, record its length
               if (val[j*3+2] == i)
               {
                  int cyc_len = pos - val[j*3+1] + 1;
                  cycdata << ',' << cyc_len;
               }
            }
            cycdata << endl;
         }

   }
   cycdata.close();
}

// Count_coprime_cycle but with a = 0
// Records the Cycle_data for trinomials in the form 
// f(x)=a+x+cx^d for every degree d coprime to p-1 over F_p, n <= p < N
// csv line format: prime, degree, 0, c, num_roots, num_periodic_pts, min_cycle_len, max_cycle_len, num_cycles\n
void Count_a0_cycle(int n, int N)
{
   vector<int> prime_list = Primes(n, N);
   fstream cycdata;
   cycdata.open("a0_cyc_data.csv",fstream::out);
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      for (int d = (int)pow(prime_list[i], 0.6); d < prime_list[i] - 1; d++)
      {
         if (!Coprime(d, prime_list[i]-1))
            continue;

         for (int c = 1; c < prime_list[i]; c++)
         {
            Cycle_data tri_data = Trinomial_cycle(prime_list[i], 0, c, d);
            cycdata << prime_list[i] << ',' << d << ',' << 0 << ',' << c << ',' << tri_data.roots << ',' << tri_data.cycle_len << ',' << tri_data.cycle_min << ',' << tri_data.cycle_max << ',' << tri_data.cycle_num << '\n';
         }

      }
      cycdata.flush();
   }  
   cycdata.close();
}

// for p between n and N
// Records the size of value set of f(x)=cx^d+x for all d >= pow(prime_list[i] coprime to p-1
// csv line format: p, d, c, value_set_size\n
void Count_a0_valueset(int n, int N)
{
   vector<int> prime_list = Primes(n, N);
   fstream cycdata;
   cycdata.open("a0_valueset.csv",fstream::out);
   for (size_t i = 0; i < prime_list.size(); i++)
   {
      for (int d = (int)pow(prime_list[i], 0.6); d < prime_list[i] - 1; d++)
      {
         if (!Coprime(d, prime_list[i]-1))
            continue;

         bool *valueset = new bool[prime_list[i]];
         for (int c = 1; c < prime_list[i]; c++)
         {
            for(int k = 0; k < prime_list[i]; k++)
               valueset[k] = false;

            int count_valueset = 0;
            for(int j = 0; j < prime_list[i]; j++)
            {
               int eval = (j + c * Pow(j, d, prime_list[i])) % prime_list[i];
               if (!valueset[eval])
               {
                  count_valueset++;
                  valueset[eval] = true;
               }
            }
            cycdata << prime_list[i] << ',' << d << ',' << c << ',' << count_valueset << '\n';
         }
         delete[] valueset;

      }
      cycdata.flush();
   }  
   cycdata.close();
}

// Prints the cycle decomposition of polynomial f(x)=a+x+cx^d explicitly 
// to standard output
// including numbers leading up to a cycle
// output line format for line with cycle: x f(x) f^2(x) ... f^m(x) f^{m-n+1}(x) || Cycle length: n
// output line format for line without cycle (when it reaches an element already appeared): x f(x) f^2(x) ... f^m(x) f^{m-n+1}(x) || 
// output last line: Number of periodic points: num_periodic_pts
void Step_trinomial(int p, int a, int c, int d)
{
   int cycle_length_sum = 0;

   // val is a p*3 2d array
   // val[i][j] is then rewritten as val[i*3+j]
   int *val = new int[p*3];
   // val[][0] being -1 means not visited
   // val[][1] is the position in cycle
   // val[][2] is i
   for (int i = 0; i < p*3; i++)
      val[i] = -1;
   
   for (int i = 0; i < p; i++)
   {
      if (val[i*3+0] != -1)
         continue;

      int j = i;
      int pos = 0;
      while (val[j*3+0] == -1)
      {
         cout << j << " ";
         pos++;
         val[j*3+0] = (a + j + c * Pow(j, d, p)) % p;
         val[j*3+1] = pos;
         val[j*3+2] = i;
         j = val[j*3+0];
      }
         
      cout << j << " || \n";

      // if the cycle is detected in this iteration, record its length
      if (val[j*3+2] == i)
      {
         int cyc_len = pos - val[j*3+1] + 1;
         cycle_length_sum += cyc_len;
         cout << "Cycle length: " << cyc_len << endl;
      }
   }
   delete [] val;
   cout << endl << "Number of periodic points: " << cycle_length_sum << endl;
}

int main(int argc, char **argv)
{   
   if ( argc == 2 )
   {
      char* p_end;  
      const long N = std::strtol( argv[ 1 ], &p_end, 10 );  
      if ( argv[ 1 ] == p_end )
      {
         fprintf ( stderr, "Usage for finding average behavior for p up to N:  <executable> N\n" );  
         return -1;  
      }
      Count_cycle_average(0, N);
   }
   else if ( argc == 5 )
   {
      int a[ 4 ];  
      for ( int i = 1;  i < argc;  i++ )
      {
         char* p_end;  
         const long element = std::strtol( argv[ i ], &p_end, 10 );  
         if ( argv[ i ] == p_end )
         {
            fprintf ( stderr, "Usage for cycle decomposition for 1 polynomial:  <executable> p(prime) a(constant) c(coefficient for x^d) d(degree)\n" );  
            return -1;  
         }
         a [ i - 1 ] = element;  
      }
      Step_trinomial(a[0], a[1], a[2], a[3]);
   }
   else if ( argc == 3 )
   {
      char* p_end;  
      const long Min = std::strtol( argv[ 1 ], &p_end, 10 );  
      if ( *argv[ 1 ] == 'p' )
      {
         const long P = std::strtol( argv[ 2 ], &p_end, 10 );  
         if ( argv[ 2 ] == p_end )
         {
            fprintf ( stderr, "Usage for counting cycle lengths for cx^d+x+{0,1} in Fp:  <executable> 'p' p\n" );  
            return -1;  
         }
         Count_large_prime_cycle(P);
         return 0;
      }
      
      if ( argv[ 1 ] == p_end )
      {
         fprintf ( stderr, "Usage for counting polynomials degree coprime to p-1 for p in [Min, N):  <executable> Min N\n" );  
         return -1;  
      }
      const long N = std::strtol( argv[ 2 ], &p_end, 10 );  
      if ( argv[ 2 ] == p_end )
      {
         fprintf ( stderr, "Usage for counting polynomials degree coprime to p-1 for p in [Min, N):  <executable> Min N\n" );  
         return -1;  
      }
      Count_coprime_cycle(Min, N);
   }
   else if ( argc == 4 )
   {
      char* p_end;  
      const long Min = std::strtol( argv[ 2 ], &p_end, 10 );  
      if ( argv[ 2 ] == p_end )
      {
         fprintf ( stderr, "Usage for counting value set of x+cx^d polynomials degree coprime to p-1 for p in [Min, N):  <executable> v Min N\n" );  
         return -1;  
      }
      const long N = std::strtol( argv[ 3 ], &p_end, 10 );  
      if ( argv[ 3 ] == p_end )
      {
         fprintf ( stderr, "Usage for counting value set of x+cx^d polynomials degree coprime to p-1 for p in [Min, N):  <executable> 'v' Min N\n" );  
         return -1;  
      }
      Count_a0_valueset(Min, N);
   }
   
   else
   {
      fprintf ( stderr, "Usage for finding average behavior for p up to N:  <executable> N\nUsage for cycle decomposition for 1 polynomial:  <executable> p(prime) a(constant) c(coefficient for x^d) d(degree)\nUsage for counting polynomials degree coprime to p-1 for p in [Min, N):  <executable> Min N\nUsage for counting value set of x+cx^d polynomials degree coprime to p-1 for p in [Min, N):  <executable> v Min N\nUsage for counting cycle lengths for cx^d+x+{0,1} in Fp:  <executable> 'p' p\n" );  
      return -1;  
   }
}