#include <math.h>
#include <vector>

using std::vector;

inline int gcd(int a, int b)
{
   while(b)
   {
      int t = a % b;
      a = b;
      b = t;
   }
   return a;
}

inline bool Coprime(int a, int b)
{
   if(!((a | b) & 1))
      return false; // Both are even numbers, divisible by at least 2.
   return 1 == gcd(a, b);
}

// computes (a * b) % mod
inline int Mult(int a, int b, int mod)
{
   int res = 0; 
   a = a % mod;
   if (a == 0) return 0; 

   while (b > 0)
   {
      if (b % 2 == 1)
         res = (res + a) % mod;

      a = (a * 2) % mod;
      b /= 2;
   }
   
   return res;
}

// Integer powers (x^exp) mod p for p up to INT_MAX/2, above that it might overflow
inline int Pow(int x, int exp, int p)
{
   int res = 1;
   x = x % p; 
   if (x == 0) return 0; 

   while (exp > 0)
   {
      if (exp % 2 == 1)
         res = Mult(res, x, p);

      x = Mult(x, x, p);
      exp /= 2; 
   }

   return res;
}

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


int primitive_root(int p, int largest_p)
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