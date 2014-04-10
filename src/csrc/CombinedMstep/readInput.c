#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include "nrutil.h"

#define LINELENGTH 1000
#define DIGITS 20

int *ReadInputMatI(FILE *file, int numRows, int numCols) {
	char buffer[LINELENGTH];

	int **matrix = (int **) imatrix(1, numRows, 1, numCols);
	int i, j, k, l;
	for (i=1; 1; i++) {
		if (fgets(buffer, LINELENGTH-1, file) == NULL)
			break;
		k = 0;
		l = 0;
		j = 1;
		while (buffer[k] != '\n') {
			while (buffer[k] != ',')
				k++;
			buffer[k] = '\0';
			matrix[i][j] = atoi(&buffer[l]);
			l = k;
			j++;

		}
		matrix[i][j+1] = atoi(&buffer[l]);

	}

	return matrix;
}

int *ReadInputI(FILE *file, int *T) {
	int lines_allocated = 128;
	char buffer[LINELENGTH];

	int *sequence = (int *) malloc(sizeof(int) * (lines_allocated + 1));
	if (sequence==NULL) {
		fprintf(stderr, "Out of memory. \n");
		exit(1);
	}
	if (file == NULL) {
		fprintf(stderr, "Error opening file. \n");
		exit(2);
	}

	int i, j;
	for (i=1; 1; i++) {
		if (i >= lines_allocated) {
			int new_size;
			new_size = lines_allocated << 1;
			sequence = (int *) realloc(sequence, sizeof(double) * (new_size + 1));
			if (sequence == NULL) {
				fprintf(stderr, "Out of memory.\n");
				exit(3);
			}
			lines_allocated = new_size;
		}
		if (fgets(buffer, LINELENGTH-1, file) == NULL)
			break;

		sequence[i] = atoi(buffer);

	}
	*T = i-1;
	sequence = (int *) realloc(sequence, sizeof(double) * (i));
	return sequence;
}

double *ReadInputD(FILE *file, int *T) {
	int lines_allocated = 128;
	char buffer[LINELENGTH];

	double *sequence = (double *) malloc(sizeof(double) * (lines_allocated + 1));
	if (sequence==NULL) {
		fprintf(stderr, "Out of memory. \n");
		exit(1);
	}
	if (file == NULL) {
		fprintf(stderr, "Error opening file. \n");
		exit(2);
	}

	int i, j;
	for (i=1; 1; i++) {
		if (i >= lines_allocated) {
			int new_size;
			new_size = lines_allocated << 2;
			sequence = (double *) realloc(sequence, sizeof(double) * (new_size + 1));
			if (sequence == NULL) {
				fprintf(stderr, "Out of memory.\n");
				exit(3);
			}
			lines_allocated = new_size;
		}
		if (fgets(buffer, LINELENGTH-1, file) == NULL)
			break;

		for (j = strlen(buffer) -1; j >= 0 && (buffer[j] == '\n' || buffer[j] == '\r'); j--);
		buffer[j+1] = '\0';
		sequence[i] = atof(buffer);

	}
	*T = i-1;
	sequence = (double *) realloc(sequence, sizeof(double) * (i));
	return sequence;
}
