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

// Integer powers (x^exp) mod p for p up to 46340, above that it might overflow
int Pow(int x, int exp, int p)
{
  if (exp == 0) return 1;
  if (exp == 1) return x;
  
  int tmp = Pow(x, exp/2, p);
  if (exp % 2 == 0) return (tmp * tmp) % p;
  else return (((x * tmp) % p) * tmp) % p;
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