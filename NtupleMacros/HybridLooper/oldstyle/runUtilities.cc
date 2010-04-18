#include "runUtilities.h"

#include <fstream>
#include "TString.h"

int highestGoodrun () {

     ifstream stream("goodruns.txt");
     TString line = "";
     Int_t runMax = 0;
     Int_t run = 0;

     while (line.ReadLine(stream)) {
        TString srun = "";
        for (Int_t i = 0; i < line.Length(); ++i) {      
                TString character = line[i];
                if (character.IsWhitespace()) break;
                srun += character;
        }
        run = srun.Atoi();           
        if (run > runMax) runMax = run;
     }

     return runMax;

}

