* for build
- gcc -fopenmp main.c bitmap.c

* for exec
- a.out [# of threads] [input image] [output image]

- image should be 24bit bmp
- With the sample main.c, [# of threads] is meaningless because the program is not parallelized.

