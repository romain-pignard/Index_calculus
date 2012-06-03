/*********

ECE 6280 project
* Romain Pignard
* Solving DLP using index calculus
* usage : ./index_calculus number 
* you obtain the 'a' such that g^a = number mod p



********/
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"



mpz_t* prime_gen(int num_prime)
{
	//generation of an array with num_prime consecutives primes 
	int i = 1;	
	mpz_t * prime_array = (mpz_t *) malloc(num_prime * sizeof(mpz_t));
	//prime_array[0..num_prime-1]
	mpz_init(prime_array[0]);
	mpz_set_str(prime_array[0],"2",10);
	for(i =1; i< num_prime; i++)
 	{	
	mpz_init(prime_array[i]);
	//nice function provided by GMP
	mpz_nextprime(prime_array[i], prime_array[i-1]);

	}
	//re-set of the first element due to unknown modification 
	mpz_set_str(prime_array[0],"2",10);
	return prime_array;
}



long int* index_array(mpz_t g_i,mpz_t exponent, mpz_t * prime_array, int num_prime)
{
	//factorization with the primes
	
	//we try to divide the number we want to factorize by the different primes.
	// we put the corresponding exponents in an array.  
	//including the "i" of g^i at the end of the row to obtain a system of linear 
	//equations
	
	
	
	
	
	long int* index_tab = (long int*) malloc((num_prime + 1) * sizeof(long int));
	//index_tab[0..num_prime]
	
	int j = 0;
	mpz_t remainder;
	mpz_init(remainder);
	mpz_t quotient;
	mpz_init(quotient);
	index_tab[num_prime] = mpz_get_si(exponent); //"i" of g^i 
	for(j = 0; j < num_prime; j++)
	{
		
	
		index_tab[j] = 0;	
		//division of the number by the prime
		mpz_cdiv_qr(quotient,remainder, g_i, prime_array[j]);
		
		 
		while((mpz_cmp_d(remainder, 0) == 0))
		{
			//while the number is divisble by the prime, 
			//we increment its exponent  
			index_tab[j]++;
			mpz_set(g_i, quotient);
			//division of the number by the prime.
			mpz_cdiv_qr(quotient,remainder, g_i, prime_array[j]);
			
		}
		
	}
	if (mpz_cmp_d(quotient, 1) != 0)
	{
		
		//we return a NULL pointer if the fatorization failed
		return NULL;
	}
	else
	{
		
		return index_tab;
	}	 
}

long int ** get_relations(int max_rel, mpz_t* prime_array, int num_prime,mpz_t g, mpz_t mod_base, mpz_t max_exp)
{
	
	
	int num_rel = 0; //number of succesful factorizations
	long int* temp_array;  
	mpz_t g_k;
	mpz_t k_exp;
	
	mpz_init(g_k);

	mpz_init(k_exp);
	
	
	long int ** relation_array = (long int **) malloc(max_rel*sizeof( long int*));
	
	
	//PRNG setup /!\ no entropy => same primes for every run
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	
	int k;
	for(k = 1; num_rel < max_rel+1; k++)
	{
		//get a new random exponent < max_exp (max_exp = q in our case)
		mpz_urandomm (k_exp,state,max_exp);
		
		mpz_powm(g_k, g, k_exp, mod_base); 
		
				
		temp_array = index_array(g_k, k_exp, prime_array, num_prime);
		if (temp_array != NULL)
		{
			//if the factorization is successful, we keep the relation
			relation_array[num_rel] = temp_array;
			num_rel++;
		}	
	} 
	return relation_array;
}
 
  void print_row(mpz_t* row,int size_row)
 {
	 int i;
	 for(i=0; i<size_row; i++)
	 {gmp_printf(" %Zd", row[i]);}
	 printf("\n");
 }
 
 
 void sub_row(mpz_t* row1, mpz_t* row2, mpz_t fact, int size_row,mpz_t mod_base)
 {
		//row1 = row1 - fact*row2
		mpz_t in_fact;
		mpz_init(in_fact);
		mpz_set(in_fact,fact);
		int i;
		mpz_t temp;
		mpz_init(temp);
		for(i=0;i<size_row;i++)
		{ 
			
			mpz_mul(temp,row2[i],in_fact);
			mpz_sub(row1[i],row1[i],temp);
			mpz_mod(row1[i],row1[i],mod_base);
			
		} 
	
	 
 }  
 
 void normalize_row(mpz_t* row, mpz_t fact_norm, int size_row,mpz_t mod_base)
 {
	int i;
	mpz_t inv_fact; 
	mpz_init(inv_fact);
	
	//inv_fact = inverse of fact_norm in the group Z[mod_base]*
	mpz_invert(inv_fact,fact_norm,mod_base);
	
	mpz_t temp; 
	mpz_init(temp);
	
	for(i=0;i<size_row;i++)
	{
		mpz_mul(temp, row[i], inv_fact); 
		mpz_mod(row[i], temp, mod_base);
	} 
 }  
 
 
 void print_array(mpz_t ** array, int size_row, int size_column)
 {
	 int i;
	 for(i=0; i<size_column; i++)
	 {print_row(array[i],size_row);}
	 printf("\n");
 }
 
 

 
 
 long int* solve(long int ** rel_array, mpz_t mod_base,int n)
 {
	//solve the linear system of equations
	mpz_t ** mp_array = (mpz_t **) malloc((n)*sizeof(mpz_t *));
	if (mp_array == NULL){printf("fail\n");}
	mpz_t * temp_line;
	int k;
	int i=0;
	int j;
	for(i=0;i<n;i++)
	{
		mp_array[i] = (mpz_t *) malloc((n+1)*sizeof(mpz_t));
		for(j=0;j<=n;j++)
		{
			mpz_init_set_ui(mp_array[i][j], rel_array[i][j]);
		}
	}	
	for(k=0;k<n;k++)
	{
		i=k;	
		while(((i<n) && mpz_cmp_d(mp_array[i][k], 0) == 0)  )
		{			
			
			if  ((i == n))
			{
					
				return NULL;
			}
			i++;		
		}	
		
		temp_line = mp_array[k];
		mp_array[k] = mp_array[i];
		mp_array[i] = temp_line;
				
		normalize_row(mp_array[k],mp_array[k][k],n+1,mod_base);
		
		for(i=0;i<n;i++)
		{
			if (i != k)
			{
				
				sub_row(mp_array[i],mp_array[k],mp_array[i][k],n+1,mod_base);
			}	
		}
		
		
		
	} 
	
	long int* result = (long int*) malloc(n*sizeof(long int));
	for(i=0;i<n-1;i++)
	{
		result[i] = mpz_get_si(mp_array[i][n]);
		
	}
	return result;
 }
 
 
 
void display_relations(int num_rel,int num_prime,long int** rel_array,mpz_t * prime_array, mpz_t g, mpz_t p)
{
		int i;
		int j;
		for(i=0; i<num_rel; i++)
		{
			gmp_printf("%Zd^%d mod %Zd = ",g,rel_array[i][num_prime],p); 
			for(j=0;j< num_prime; j++)
			{
				if (rel_array[i][j] != 0)
				{
					gmp_printf("(%Zd^%d)*",prime_array[j],rel_array[i][j]);
				}	
			} 
			printf("1\n");
		}
	
	
	
} 





int find_exponent(long int * log_primes,mpz_t * tab,  int num_primes, mpz_t g, mpz_t g_alpha, mpz_t q, mpz_t p)
{
	
	//final function that does 
	mpz_t g_temp_s;
	mpz_init(g_temp_s);
	
	
	mpz_t exponent;
	mpz_init(exponent);
	
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	
	mpz_t result;
	mpz_init(result);
	
	mpz_t fact;
	mpz_t log_prime; 
	mpz_init(fact);
	mpz_init(log_prime);
	
	long int * factorization = NULL; //no factorization yet
	
	while(factorization == NULL)
	{
			mpz_urandomm (exponent,state,q);
			mpz_powm(g_temp_s,g,exponent,p);
			
			mpz_mul(g_temp_s,g_temp_s,g_alpha);
			mpz_mod(g_temp_s,g_temp_s,p);
			
			factorization = index_array(g_temp_s,exponent, tab, num_primes);
			
	}
	
	int i;
	for(i=0; i< num_primes;i++)
	{
		//result = result + factorization[i] * log_primes[i]
		mpz_set_si(fact,  factorization[i]);
		mpz_set_si(log_prime,log_primes[i]);
		
		mpz_addmul(result,fact,log_prime);
		mpz_mod(result,result,q);
	}	
	return mpz_get_si(result) - mpz_get_si(exponent) ;
	
	
}

int main(int argc, char ** argv)
{
 mpz_t p;
 mpz_t q;
 mpz_t k;
 mpz_t g;
 mpz_t g_alpha;
 mpz_t l;	

 
 mpz_init(p);
 mpz_init(q);
 mpz_init(k);
 mpz_init(g);
 mpz_init(g_alpha);
 mpz_init(l);

 
 mpz_set_str(p, "2382933803", 10);
 mpz_set_str(q, "10930889", 10);
 mpz_set_str(k, "1", 10);
 mpz_set_str(g, "2084483647", 10);
 mpz_set_str(g_alpha, argv[1],10);
 mpz_set_str(l, "15", 10);
 
 int nb_rel = 11; //we choose 11 for speed.
 mpz_t *  tab = prime_gen(nb_rel);
 long int ** rel_array = get_relations(nb_rel, tab, nb_rel, g, p,q);
 
 printf("\n\nwe use these relations to obtain the result\n\n");
 display_relations(nb_rel, nb_rel,rel_array,tab,g,p); 
 printf("\n");
long int * log_primes =  solve(rel_array,q,nb_rel);

long int log = find_exponent(log_primes, tab,nb_rel, g,g_alpha,q,p);  
printf("\n\nresult = %ld\n\n", log); 
	
 
 
 /* free used memory (or not) */
 mpz_clears(p, q, k, g, NULL);
 return EXIT_SUCCESS;
}
