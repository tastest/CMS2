double Q(double A, double B, double C){
  double q = A/C*((1-2*A/B)**2);
  printf("P: %0.1f%%, Eff: %0.1f%%, Q: %0.1f%%\n", A/B*100, A/C*100, q*100);
  return q;
}
