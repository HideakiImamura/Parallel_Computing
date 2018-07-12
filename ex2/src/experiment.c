#include "bitmap.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[])
{
  if(argc != 2){
    fprintf(stderr, "Usage: program <# of image>");
    exit(1);
  }

  Image *colorimg;
  int m_nthread = omp_get_num_procs();
  char buf[6];
  int n, nthread, iter_n = 5;
  FILE *fp;
  
  printf("# of processor: %d\n", m_nthread);
  fp = fopen("result.csv", "a");
  fprintf(fp, "nthread,height,width,time\n");

  for(n=0; n<atoi(argv[1]); n++){
    snprintf(buf, 6, "%d.bmp", n);
    if((colorimg = Read_Bmp(buf)) == NULL){
      exit(1);
    }

    for(nthread=1; nthread<=m_nthread; nthread++){
      int i, t, maxp = colorimg->height * colorimg->width;
      double st, en, sum = 0.0;
      omp_set_num_threads(nthread);
      for(t=0; t<iter_n; t++){
        st = omp_get_wtime();
#pragma omp parallel for
        for(i=0; i<maxp; i++){
          colorimg->data[i].b = 255 - colorimg->data[i].b;
          colorimg->data[i].g = 255 - colorimg->data[i].g;
          colorimg->data[i].r = 255 - colorimg->data[i].r;
        }
        en = omp_get_wtime();
        sum += (en-st)*1000.0;
      }
      sum /= iter_n;
      fprintf(fp, "%d,%d,%d,%f\n", nthread, colorimg->height, colorimg->width, sum);
    }
    Free_Image(colorimg);
  }
  fclose(fp);
  return 0;
}
