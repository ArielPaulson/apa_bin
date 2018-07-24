#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

int main(int argc, char *argv[]) { 
	
	char *ptr;
	int COLS, CMAX;
	
	COLS = strtol(argv[2], &ptr, 10);
	CMAX = strtol(argv[3], &ptr, 10);
	
	char line[5000000];
	char thisline[5000000];
	char *dtoken;
	char *token;
	
	int irow = 0;
	int icol = 0;
	int i = 0;
	
	char *transposed[COLS];
	for (i=0; i<COLS; i++) {
		transposed[i] = (char *)malloc(sizeof(char) * CMAX);
	}
	
	FILE *infile = fopen(argv[1], "r");
	
	while (fgets(line, sizeof(line), infile) != NULL) {
		
		/* Get current line.
		   This may not be the entire line, if lines are >5M chars */
		sscanf(line, "%[^\n]", thisline);
		
		icol = 0;  // in col == trans row
		token = "1";  // initialize non-NULL
		
		while ( token != NULL ) {
			if (icol == 0) {
				// start of new line
				token = strtok(thisline, "\t");
			} else {
				token = strtok(NULL, "\t");
			}
			if (irow == 0) {
				// very first row; initializes all output rows
				strcpy(dtoken, token);
			} else {
				strcpy(dtoken, "\t");
				strcat(dtoken, token);
			}
			//printf("IROW: %d | ICOL: %d\n", irow, icol);
			strcat(transposed[icol], dtoken);  // grow output row one field at a time
			printf(" TOKEN 1: '%s' | DTOKEN 1: '%s'\n", token, dtoken);
			icol++;
		}
		
		puts("Survived!\n");
		
		/* Test if 'line' ended with a newline.
		   If so, then increment input row count */
		regex_t regex;
		int reti;
		reti = regcomp(&regex, "\n$", 0);
		if (reti) {
			fprintf(stderr, "Could not compile regex\n");
			exit(1);
		}
		reti = regexec(&regex, line, 0, NULL, 0);
		if (!reti) {
			irow++;  // EOL found; increment
		}
	}
	
	for (i=0; i<=icol; i++) {
		printf("%s\n", transposed[i]);
	}
	
	return(0);
	
}
