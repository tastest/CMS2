#include "tauify.C"

void test(){
  Tauify *tau = new Tauify("decay.txt", false);

  // test object access by index (pid, P, costheta)
  for(unsigned int i=0; i<tau->TauSize(); i++){
    cout << i << "\t" << tau->First(i) <<  "\t" << tau->Second(i) << "\t" << tau->Third(i) << endl;
    if( i > 9) exit(1);
  }

}
