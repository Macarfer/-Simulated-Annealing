#include "fileManager.h"

void openFile(const char * path,int matrix[DIMENSION][DIMENSION]){
  int i=0,j=0;
  int buff;

  FILE *f=NULL;

  if (!(f=fopen(path,"r"))){
    perror("Can't open file!");
    exit(1);
  }

for (i=0;i<DIMENSION;i++) {
  for(j=0;j<i;j++){
  fscanf(f, "%d", &buff);
  matrix[i][j]=buff;
  matrix[j][i]=buff;
}
}
  fclose(f);
}

double * openRandoms(const char * path){
  int i=0,lines=0;
  double * randoms, buff;
  FILE *f=NULL;

  if (!(f=fopen(path,"r"))){
    perror("Can't open file!");
    exit(1);
  }
  for(;fscanf(f, "%lf", &buff)!=EOF;)
    lines+=1;

  fseek(f,0,SEEK_SET);
  randoms=malloc(lines*sizeof(double));

  for(i=0;fscanf(f, "%lf", &buff)!=EOF;i++){
    randoms[i]=buff;
  }
  fclose(f);
  return randoms;
}
