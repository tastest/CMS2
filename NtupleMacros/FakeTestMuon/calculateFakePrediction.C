#include <iomanip>
#include <iostream>
void calculateFakePrediction(double nFakeables, double nFakeables_error, double fakerate, double fakerate_error) {
double prediction = nFakeables * fakerate;
double prediction_error = sqrt(nFakeables*nFakeables*fakerate_error*fakerate_error+fakerate*fakerate*nFakeables_error*nFakeables_error);
cout << "Fake prediction: " << fixed << setprecision(1) << prediction << "+-" << prediction_error << endl;
}
