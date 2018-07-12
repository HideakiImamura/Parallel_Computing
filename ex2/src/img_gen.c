#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bitmap.h"
#include <math.h>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: program <# of image>\n");
    exit(1);
  }

  int i, base_height = 100, base_width = 60, height, width;
  char buf[6];
  Image *img;
  for(i = 0; i < atoi(argv[1]); i++){
    height = base_height * pow(2, i);
    width  = base_width  * pow(2, i);
    snprintf(buf, 6, "%d.bmp", i);
    img = Create_Image(height, width);
    Write_Bmp(buf, img);
  }

  return 0;
}
