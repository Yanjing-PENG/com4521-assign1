#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "NBody.h"
#include "NBodyVisualiser.h"

#define USER_NAME "acp18aaf"		//replace with your username

void print_help();
void step(void);

struct nbody* nbodies;

N = 1000;
D = 10;
const float EPSILON = 0.00001;

int main(int argc, char *argv[]) {

	int i, j;

	double begin, end;
	double seconds;

	void (*f_p)(void);
	f_p = &step;

	//TODO: Processes the command line arguments
		//argc in the count of the command arguments
		//argv is an array (of length argc) of the arguments. The first argument is always the executable name (including path)

	//TODO: Allocate any heap memory
		//allocate memory to this variables;
	nbodies = (struct nbody*)malloc(sizeof(struct nbody) * N);
	
	//TODO: Depending on program arguments, either read initial data from file or generate random data.
	begin = omp_get_wtime();

	FILE* f = NULL;
	f = fopen("example_input.nbody", "r"); //read input file
	char line[100];
	char split[] = ",";
	char* result = NULL;

	srand((unsigned)time(NULL));

	if (f == NULL) {
		fprintf(stderr, "Error: Could not find input file \n");
		exit(1);
	}
	j = 0;
	while (!feof(f))
	{
		fgets(line, 99, f);
		if (line[0] == '#') {
			continue;
		}
		result = strtok(line, split);
		i = 0;
		while (result != NULL) {
			switch (i++)
			{
			case 0:
				nbodies[j].x = atof(result);
				if (nbodies[j].x <= EPSILON) nbodies[j].x = rand() / 32767.0;
				break;
			case 1:
				nbodies[j].y = atof(result);
				if (nbodies[j].y <= EPSILON) nbodies[j].y = rand() / 32767.0;
				break;
			case 2:
				nbodies[j].vx = atof(result);
				if (nbodies[j].vx <= EPSILON) nbodies[j].vx = 0.0;
				break;
			case 3:
				nbodies[j].vy = atof(result);
				if (nbodies[j].vy <= EPSILON) nbodies[j].vy = 0.0;
				break;
			case 4:
				nbodies[j].m = atof(result);
				if (nbodies[j].m <= EPSILON) nbodies[j].m = 1.0 / N;
				j++;
				break;
			}
			result = strtok(NULL, split);
		}
	}
	if (j != N) {
		fprintf(stderr, "Error: The number of nbodies is wrong, change the value of N. \n");
		exit(1);
	}
	if (D <= 0) {
		fprintf(stderr, "Error: The grid size is wrong, change the value of D. \n");
		exit(1);
	}


	fclose(f);

	//TODO: Depending on program arguments, either configure and start the visualiser or perform a fixed number of simulation steps (then output the timing results).

	initViewer(N, D, OPENMP, f_p);

	end = omp_get_wtime();
	seconds = (end - begin);
	printf("Calculation complete in %.2f seconds\n", seconds);

	for (i = 0; i < N; i++) {
		setNBodyPositions(&nbodies[i]);
		startVisualisationLoop();
	}


	destroyViewer();

	//free memory;
	free(nbodies);

	return 0;
}

void step(void)
{
	//TODO: Perform the main simulation of the NBody system

	float* vx, * vy, * fx, * fy, * ax, * ay, * num;
	int i, j;
	vx = (float*)malloc(sizeof(float) * N);
	vy = (float*)malloc(sizeof(float) * N);
	fx = (float*)malloc(sizeof(float) * N);
	fy = (float*)malloc(sizeof(float) * N);
	ax = (float*)malloc(sizeof(float) * N);
	ay = (float*)malloc(sizeof(float) * N);
	num = (float*)malloc(sizeof(float) * D * D);
	memset(fx, 0, sizeof(float) * N);
	memset(fy, 0, sizeof(float) * N);
	memset(num, 0, sizeof(float) * D * D);


#pragma omp parallel for default(none) shared(nbodies) private(j, i)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			vx[i] = sqrt(pow(nbodies[i].x - nbodies[j].x, 2) + pow(nbodies[j].x - nbodies[i].x, 2)); //calculate the speed in x axis
			vy[i] = sqrt(pow(nbodies[i].y - nbodies[j].y, 2) + pow(nbodies[j].y - nbodies[i].y, 2)); //calculate the speed in y axis
			fx[i] += (nbodies[j].m * (nbodies[j].x - nbodies[i].x)) / pow((pow(vx[i], 2) + pow(SOFTENING, 2)), 3 / 2);
			fy[i] += (nbodies[j].m * (nbodies[j].y - nbodies[i].y)) / pow((pow(vy[i], 2) + pow(SOFTENING, 2)), 3 / 2);

		}
		fx[i] *= G * nbodies[i].m; //calculate force in x axis
		fy[i] *= G * nbodies[i].m; //calculate force in y axis
		ax[i] = fx[i] / nbodies[i].m; //calculate acceleration in x axis
		ay[i] = fy[i] / nbodies[i].m; //calculate acceleration in y axis

		nbodies[i].x = nbodies[i].x + dt * nbodies[i].vx;
		nbodies[i].y = nbodies[i].y + dt * nbodies[i].vy;
		nbodies[i].vx = nbodies[i].vx + dt * ax[i];
		nbodies[i].vy = nbodies[i].vy + dt * ay[i];
		if (nbodies[i].x >= 0 && nbodies[i].x < D && nbodies[i].y >= 0 && nbodies[i].y < D) { // consider the border
			int location = (int)(floor(nbodies[i].x * D) + floor(nbodies[i].y * D) * D);
			if (location >= 0 && location < D * D) { //make sure the location exist
				num[location] += 1; //counter for nbodies
			}
		}
	}

	for (i = 0; i < D * D; i++) {
		num[i] = num[i] * D / N;
		if(num[i + 1] < EPSILON && i + 1 < D * D) i++;
	}

	setHistogramData(num);

	free(vx);
	free(vy);
	free(fx);
	free(fy);
	free(ax);
	free(ay);
	free(num);
}



void print_help(){
	printf("nbody_%s N D M [-i I] [-i input_file]\n", USER_NAME);

	printf("where:\n");
	printf("\tN                Is the number of bodies to simulate.\n");
	printf("\tD                Is the integer dimension of the activity grid. The Grid has D*D locations.\n");
	printf("\tM                Is the operation mode, either  'CPU' or 'OPENMP'\n");
	printf("\t[-i I]           Optionally specifies the number of simulation iterations 'I' to perform. Specifying no value will use visualisation mode. \n");
	printf("\t[-f input_file]  Optionally specifies an input file with an initial N bodies of data. If not specified random data will be created.\n");
}
