
#include "EffH1F.h"
#include <string>
#include <math.h>

class EffMulti {

	public:
		EffMulti() {}
		EffMulti(bool lessThan, float thresholdEB, float thresholdEE,
			 std::string source, std::string var, std::string det);
		~EffMulti();

		void Fill(float value, float pt, float eta, float phi, float weight);

	private:
		EffH1F *e1_pt_;
		EffH1F *e1_eta_;
		EffH1F *e1_phi_;

		bool lessThan_;
		float thresholdEB_;
		float thresholdEE_;


};

