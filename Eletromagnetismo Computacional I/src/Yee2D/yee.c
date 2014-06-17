/*
 * The MIT License (MIT)
 * Copyright (c) 2014 Leonardo Koch Kewitz
 */
#include <stdio.h>
#include <string.h>
/* #define I(x,y,cols) ((y*cols)+x) */

int I (int x, int y, int cols) {
	int index = (y*cols)+x;
	return index;
};

int calcH (int k, int skip, unsigned int sx, unsigned int sy, double CH, double * Ez, double * Hx, double * Hy, double * oldHx, double * oldHy) {
	int t,x,y,in;
	for(y = 0; y < sy-1; y++)
		for(x = 0; x < sx-1; x++)
		{
			in = I(x,y,sx);
			Hy[in] = CH * (Ez[I(x,y+1,sx)] - Ez[I(x,y,sx)]);
			Hx[in] = -CH * (Ez[I(x+1,y,sx)] - Ez[I(x,y,sx)]);
			if(k > 0)
			{
				Hy[in] += oldHy[in];
				Hx[in] += oldHx[in];
			}
		}

	return 1;
};

int calcE (int skip, int sx, int sy, double CEx, double CEy, double * Ez, double * Hx, double * Hy, double * oldEz) {
	int t,x,y,in;
	for(y = 1; y < sy-1; y++)
		for(x = 1; x < sx-1; x++)
		{
			in = I(x,y,sx);
			Ez[in] = oldEz[in] + (CEx * ( Hy[I(x,y,sx)] - Hy[I(x,y-1,sx)] )) - (CEy * (Hx[I(x,y,sx)] - Hx[I(x-1,y,sx)]));
		}

	return 1;
};
