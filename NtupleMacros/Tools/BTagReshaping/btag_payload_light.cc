#include "btag_payload_light.h"

float mistag_CSVL(float eta,float x, float scale) {
  float mean=-1.,min=-1.,max=-1.;
  if(eta > 0.0 && eta <= 0.5) //   if( Atagger == "CSVL" && sEtamin == "0.0" && sEtamax == "0.5")
      { 
        mean = ((1.04901+(0.00152181*x))+(-3.43568e-06*(x*x)))+(2.17219e-09*(x*(x*x)));
        min = ((0.973773+(0.00103049*x))+(-2.2277e-06*(x*x)))+(1.37208e-09*(x*(x*x)));
        max = ((1.12424+(0.00201136*x))+(-4.64021e-06*(x*x)))+(2.97219e-09*(x*(x*x)));
      } 
  if(eta > 0.5 && eta <= 1.0) //    if( Atagger == "CSVL" && sEtamin == "0.5" && sEtamax == "1.0")
      {
        mean = ((0.991915+(0.00172552*x))+(-3.92652e-06*(x*x)))+(2.56816e-09*(x*(x*x)));
        min = ((0.921518+(0.00129098*x))+(-2.86488e-06*(x*x)))+(1.86022e-09*(x*(x*x)));
        max = ((1.06231+(0.00215815*x))+(-4.9844e-06*(x*x)))+(3.27623e-09*(x*(x*x)));
      } 
  if(eta > 1.0 && eta <= 1.5) //    if( Atagger == "CSVL" && sEtamin == "1.0" && sEtamax == "1.5")
      {
        mean = ((0.962127+(0.00192796*x))+(-4.53385e-06*(x*x)))+(3.0605e-09*(x*(x*x)));
        min = ((0.895419+(0.00153387*x))+(-3.48409e-06*(x*x)))+(2.30899e-09*(x*(x*x)));
        max = ((1.02883+(0.00231985*x))+(-5.57924e-06*(x*x)))+(3.81235e-09*(x*(x*x)));
      }
  if(eta > 1.5 && eta <= 2.4) //    if( Atagger == "CSVL" && sEtamin == "1.5" && sEtamax == "2.4")
      {
        mean = ((1.06121+(0.000332747*x))+(-8.81201e-07*(x*x)))+(7.43896e-10*(x*(x*x)));
        min = ((0.983607+(0.000196747*x))+(-3.98327e-07*(x*x)))+(2.95764e-10*(x*(x*x)));
        max = ((1.1388+(0.000468418*x))+(-1.36341e-06*(x*x)))+(1.19256e-09*(x*(x*x)));
      }
 if(scale > 0)     return (mean + (max-mean)*scale);
  if(scale< 0)     return (mean + (mean-min)*scale);
   return mean; 
}

float mistag_CSVM(float eta,float x, float scale) {
  float mean=-1.,min=-1.,max=-1.;
  if(eta > 0.0 && eta <= 0.8) //  if( Atagger == "CSVM" && sEtamin == "0.0" && sEtamax == "0.8")
      { 
        mean = ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)));
        min = ((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)));
        max = ((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)));
      }
  if(eta > 0.8 && eta <= 1.6) //    if( Atagger == "CSVM" && sEtamin == "0.8" && sEtamax == "1.6")
      {
        mean = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
        min = ((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)));
        max = ((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)));
      } 
  if(eta > 1.6 && eta <= 2.4) //    if( Atagger == "CSVM" && sEtamin == "1.6" && sEtamax == "2.4")
      {
        mean = ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)));
        min = ((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)));
        max = ((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)));
      }
 if(scale > 0)     return (mean + (max-mean)*scale);
  if(scale< 0)     return (mean + (mean-min)*scale);
   return mean; 
}

float mistag_CSVT(float eta,float x, float scale) {
  float mean=-1.,min=-1.,max=-1.;
  if(eta > 0.0 && eta <= 2.4) //    if( Atagger == "CSVT" && sEtamin == "0.0" && sEtamax == "2.4")
      {
        mean = ((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x)));
        min = ((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)));
        max = ((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x)));
      }
  if(scale > 0)     return (mean + (max-mean)*scale);
  if(scale< 0)     return (mean + (mean-min)*scale);
   return mean; 
}
