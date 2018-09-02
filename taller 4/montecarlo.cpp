/* C++ program for estimation of Pi using Monte
   Carlo Simulation */
#include <string>
#include <math.h>
#include <time.h>
// Defines precision for x and y values. More the
// interval, more the number of significant digits
// #define interval 10000
using namespace std;

double calcPiMontecarlo(long interval){
	double rand_x, rand_y, origin_dist, pi;
	long circle_points = 0;
	// Total Random numbers generated = possible x
	// values * possible y values
	for (long i = 0; i < interval; i++) {
 
		// Randomly generated x and y values
		rand_x = double(rand() % (interval + 1)) / interval;
		rand_y = double(rand() % (interval + 1)) / interval;
 
		// Distance between (x, y) from the origin
		origin_dist = rand_x * rand_x + rand_y * rand_y;
 
		// Checking if (x, y) lies inside the define
		// circle with R=1
		if (origin_dist <= 1.0){
			circle_points++;
		}
			
 

 
		// estimated pi after this iteration
		
 
		// For visual understanding (Optional)
		//cout << rand_x << " " << rand_y << " " << circle_points
		//     << " " << square_points << " - " << pi << endl << endl;
 
		// Pausing estimation for first 10 values (Optional)
		/*if (i < 20)
			getchar();*/
	}
	return double(4 * circle_points) / interval;


}


int main(int argc, char const *argv[]){
	if(argc != 2){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
	long interval=stol(argv[1], nullptr);
	double pi;
	
	// Initializing rand()
	srand(time(NULL));
	pi=calcPiMontecarlo(interval);
	
 
	// Final Estimated Value
	//cout << "\nFinal Estimation of Pi = " << pi;
	printf("Approximation of pi after %ld tosses %.9f %.9f %% \n", interval, pi, fabs(M_PI-pi)*100/M_PI);
	return 0;
}