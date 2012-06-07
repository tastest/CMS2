#ifndef SS2012ANALYSISSKIM_H
#define SS2012ANALYSISSKIM_H

#include <string>

void ntupleFilterSameSign
(
    const std::string& infile,
    const std::string& outfile,
    const std::string& jetcorr_path,
	bool do20_20 = false, 
    bool btag2 = false,
    bool jets2 = false,
    bool printPass = false
);

#endif // SS2012ANALYSISSKIM_H
