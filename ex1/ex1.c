#include<omp.h>
#include<stdio.h>

int main() {
  char a[11]={'h','e','l','l','o',' ','w','o','r','l','d'};
  int b[11]={0,0,0,0,0,0,0,0,0,0};
  int i;

#pragma omp parallel for
  for(i=0;i<11;i++) {
    b[i]=a[i];
    printf("%c, スレッド番号=%d, スレッド数=%d\n",
              b[i], omp_get_thread_num(), omp_get_num_threads());
  }
  
  for(i=0;i<11;i++)
    printf("b[%d]=%c\n",i,b[i]);
  
  return 0;
}
