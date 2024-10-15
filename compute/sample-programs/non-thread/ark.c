
public int main() {

  public int K = 3;
  public int aprime[K]; // uniformly random vector in Z_Q
  private int sprime[K]; // secret key, bit vector, shared beforehand
  public int m = 1; // message 

  public int e = 1; // discrete gaussian noise

  smcinput(aprime, 1, K); 
  smcinput(sprime, 1, K); 

  public int i;
  for(i = 0; i < K; i++){
    printf("s_prime[%i] %d\n", i, smcopen(sprime[i]));
  }


  private int a; 
  a = aprime @ sprime; 
    printf("a_prime dot s_prime %d\n",  smcopen(a));



  return 0;
}
