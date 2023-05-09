// Sai Bhargav Mandavilli
// EEL6763

/* sed.c */
/* based on http://www.doc.gold.ac.uk/~mas02fl/MSC101/ImageProcess/edge.html */

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define WIDTH 5000
#define HEIGHT 5000

int main(int argc, char* argv[])
{
	FILE *imgInput, *imgOutput; // file pointers to store input and output image text
	unsigned int X, Y;
	int	I, J, SUM, width, height, K, L;
	long sumX, sumY;
	int *inputText = (int *)malloc(sizeof(int) * (WIDTH*HEIGHT)); // array to store pixel values from input text file
	int *edgeImage = (int *)malloc(sizeof(int) * (WIDTH*HEIGHT)); // stores chunks of image pixel values
	
	int GX[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1}; // GX sobel mask
	int GY[9] = {1, 2, 1, 0, 0, 0, -1, -2, -1}; // GY sobel mask
	
	// read pixel values from input text file
	imgInput = fopen("input.txt", "r");
	for(L=0; L<HEIGHT; L++){
		for(K=0; K<WIDTH; K++)
			fscanf(imgInput, "%d", &inputText[L*WIDTH + K]);
	}
	fclose(imgInput);
		
	for (Y = 0; Y <= HEIGHT - 1; Y++){
        for (X = 0; X <= WIDTH - 1; X++){
            sumX = 0;
            sumY = 0;
			
			/* image boundaries */
            if (Y == 0 || Y == HEIGHT - 1)
                SUM = 0;
            else if (X == 0 || X == WIDTH - 1)
                SUM = 0;
			
			/* Convolution starts here */
			else{
				/*-------X GRADIENT APPROXIMATION------*/
                for (I = 0; I < 3; I++)
                {
                    for (J = 0; J < 3; J++)
                    {
                        int t = 3 * I + J;
                        sumX = sumX + (int)((*(inputText + X + I + (Y + J) * WIDTH)) * GX[t]);
                    }
                }
				/*-------Y GRADIENT APPROXIMATION-------*/
                for (I = 0; I < 3; I++)
                {
                    for (J = 0; J < 3; J++)
                    {
                        int t = 3 * I + J;
                        sumY = sumY + (int)((*(inputText + X + I + (Y + J) * WIDTH)) * GY[t]);
                    }
                }
				/*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
				SUM = abs(sumX) + abs(sumY);
			}
			if(SUM > 255) // for 8-bit precision
				SUM = 255;
			else if(SUM < 0)
				SUM = 0;
			edgeImage[(Y*WIDTH) + X] = (unsigned int)(SUM);
		}
	}
	
	imgOutput = fopen("processed_matrix.txt", "w+"); 	// output text file
	for(L=0; L<HEIGHT; L++){
		for(K=0; K<WIDTH; K++)
			fprintf(imgOutput, "%d\t", edgeImage[L*WIDTH + K]);
		fprintf(imgOutput, "\n");
	}
	fclose(imgOutput);
	
	return 0;
}
