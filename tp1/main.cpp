#include <iostream>
#include <math.h> 

/*
 * pour la semaine prochaine
 * 
 * newton
 * 
 * p = z^3 - 1
 * 
 * dp(z) = p'
 * 
 * lim zn = zn-1 - p(zn-1) / dp(zn-1)
 * 
 * bifurcations ????
 * 
 */





using namespace std;

void derive() {
	
	int n = 100;
	double xstart = 1;
	double xend = 0;
	
	for (int i = 0; i < n; i++) {

		double t = (double)i / (double)n;
		double x = xstart + (xend - xstart) * t;

		float derive = (cos(1 + x) - cos(1)) / x + sin(1);
		float derive2 = (cos(1 + x) - cos(1 - x)) / (2*x) + sin(1);

		cout << x << " " << derive << " " << derive2 << endl;
	}
	
	
}





int main() {
	
	derive();
	
}










// cmake --build /Users/frenkield/projects/sorbonne/5MM30-edp/test1/cmake-build-debug --target test1 -- -j 4

// cmake .
// cmake --build .
