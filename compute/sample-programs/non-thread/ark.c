
public int main() {

  public int K = 1024;
  public int<128> aprime[K]; // uniformly random vector in Z_Q
  public int<128> bprime; // computed ctx
  private int<128> sprime[K]; // secret key, bit vector, shared beforehand

  smcinput(sprime, 1, K); 
  smcinput(aprime, 2, K); 
  smcinput(bprime, 2); 

  return 0;
}
