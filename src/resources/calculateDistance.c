#include <stdio.h>
#include "fileManager.h"
int main(int argc,char * argv[]){
    int matrix[100][100];
    int distance=0;
    int i=0;
    openFile("distances_100.txt",matrix);
    distance+=matrix[0][atoi(argv[1])];
    //printf("%d\t%d\n",0,atoi(argv[1]));
    for (i=1;i<(argc-1);i++){
        distance+=matrix[atoi(argv[i])][atoi(argv[i+1])];
        //printf("%d\t%d\n",atoi(argv[i]),atoi(argv[i+1]));
    }
       
        distance+=matrix[atoi(argv[99])][atoi(argv[0])];
       // printf("%d\t%d\n",atoi(argv[99]),atoi(argv[0]));
        printf("Distancia: %d\n",distance);
}