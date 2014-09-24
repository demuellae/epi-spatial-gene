#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "nrutil.h"
#include "hmmTree.h"

#define BUFFERSIZE 5000
#define NUMGENES 1   //5510
#define NUMTISSUES 811


int main() {
  char buffer[BUFFERSIZE];
  char *filename = "../../../../matrix_1224667147767.txt";
  int geneIndex = -1;
  int i;

  double *O = dvector(1, NUMTISSUES);
  double *P = ivector(1, NUMTISSUES);
  
  FILE *matrix = fopen(filename, "r");
  FILE *z = fopen("z_values", "w");

  
  while(fgets(buffer, BUFFERSIZE, matrix) != NULL) {
    if (geneIndex != -1) {
      i = 0;
      while (buffer[i] != ',') { /* skip the first col*/
	i++;
      }
      i++;
      ReadCommaSequence(&buffer[i], O);
      
    }
    if (geneIndex == NUMGENES) {
      break;
    }
    geneIndex++;
  }

  return 0;
}
