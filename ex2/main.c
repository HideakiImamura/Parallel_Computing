#include"bitmap.h"
#include<stdlib.h>
#include<stdio.h>
#include<omp.h>

int main(int argc, char *argv[])
{
	if(argc != 4){
		fprintf(stderr, "Usage: program <# of threads> <inputfile> <outputfile>\n");
		exit(1);
	}

	Image *colorimg;
	int nthread = atoi(argv[1]);

	if((colorimg = Read_Bmp(argv[2])) == NULL){
		exit(1);
	}

	// In colorimag, upperleft (0,0), 
	// lowerright (colorimg->width-1,colorimg->height-1)
	// For position (i,j), colorimg->data[j*colorimg->width+i]
	int i, maxp = colorimg->height * colorimg->width;
	double st, en;

	omp_set_num_threads(nthread);
	st = omp_get_wtime();
#pragma omp parallel for
	for(i=0; i<maxp; i++){
		colorimg->data[i].b = 255 - colorimg->data[i].b;
		colorimg->data[i].g = 255 - colorimg->data[i].g;
		colorimg->data[i].r = 255 - colorimg->data[i].r;
	}
	en = omp_get_wtime();

	if(Write_Bmp(argv[3], colorimg)){
		exit(1);
	}

	Free_Image(colorimg);

	printf("# of threads: %d\n", nthread);
	printf("# of processors: %d\n", omp_get_num_procs());
	printf("Elapsed time: %g [ms]\n", (en-st)*1000.0);

	return 0;
}

