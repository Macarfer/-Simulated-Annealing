#include <stdio.h>
#include "simulatedAnnealing.h"

int main(int argc, char const *argv[]) {

  switch (argc){
    case 2:
    initWithoutRandom(argv[1]);
    generateGreedyInitialSolution();
    run();
    break;
    case 3:

    initWithRandom(argv[1],argv[2]);
    generateGreedyInitialSolution();
    run();
    break;
    default:
      printf("Error on function usage!\n Usage:\ttabuSearch pathToDistancesFile [pathToRandomNumbers]\n");
    break;
  }
}
