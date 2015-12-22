/**
    CSC: Circular Sequence Comparison
    Copyright (C) 2015 Solon P. Pissis, Ahmad Retha, Fatima Vayani 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <float.h>
#include <sys/time.h>

#include "csc.h"
#include "sacsc.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

int EDNA[90];
int BLOSUM[91];

void init_substitution_score_tables ()
{
    int i;
    char edna[] = "ATGCSWRYKMBVHDN";
    for ( i = 0; i < 15; i ++ ) {
	EDNA[(int)edna[i]] = i;
    }
    EDNA[(int)'U'] = 1; //Setting RNA U=T
    char blosum[] = "ARNDCQEGHILKMFPSTWYVBZX*";
    for ( i = 0; i < 24; i ++ ) {
	BLOSUM[(int)blosum[i]] = i;
    }
}

double delta ( char a, char b, char * alphabet )
{
    if ( a == DEL || b == DEL ) {
	return 0;
    }
    if ( strcmp ( alphabet, ALPHABET_PROT ) == 0 )
    {
	return ( double ) EBLOSUM62_matrix[ BLOSUM[(int)a] ][ BLOSUM[(int)b] ];
    }
    else
    {
	return ( double ) EDNAFULL_matrix[ EDNA[(int)a] ][ EDNA[(int)b] ];
    }
}

unsigned int nw ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, double o, double e, double * score, char * alphabet )
{
	int i, j;
	double g = o;
	double h = e;
	double u, v, w;

	double ** T;
	double ** I;
	double ** D;

	if ( ( T = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( I = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: I could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( I[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: I could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( D = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: D could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( D[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: D could not be allocated!\n");
                        return ( 0 );
                }
        }

	init_substitution_score_tables ();

        for ( i = 0; i < m + 1; i++ )
	{	
		D[i][0] = -DBL_MAX;
		I[i][0] = -DBL_MAX;
	}

        for ( j = 1; j < n + 1; j++ )
	{	
		D[0][j] = -DBL_MAX;
		I[0][j] = -DBL_MAX;
	}

	T[0][0] = 0; 

	if ( m > 0 )
		T[1][0] = g; 

        for ( i = 2; i < m + 1; i++ )
		T[i][0] = T[i - 1][0] + h;

	if ( n > 0 )
		T[0][1] = g;

        for ( j = 2; j < n + 1; j++ )
		T[0][j] = T[0][j - 1] + h;

	for( i = 1; i < m + 1; i++ )
	{
        	for( j = 1; j < n + 1; j++ )
        	{
			D[i][j] = cscmax ( D[i - 1][j] + h, T[i - 1][j] + g );
			u = D[i][j];

			I[i][j] = cscmax ( I[i][j - 1] + h, T[i][j - 1] + g );
			v = I[i][j];

			w = T[i - 1][j - 1] + delta ( t[j - 1], p[i - 1], alphabet );

			T[i][j] = cscmax ( w, cscmax ( u, v ) );
        	}
    	}
	( * score ) = T[m][n];

        for ( i = 0; i < m + 1; i ++ )
	{
		free ( D[i] );
		free ( I[i] );
		free ( T[i] );
	}
	free ( I );
	free ( D );
	free ( T );
	
	return EXIT_SUCCESS;
}

unsigned int sacsc_refinement (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance )
{
	unsigned int rot;
	unsigned int dist;
	circular_sequence_comparison ( x, y, sw, &rot, &dist );

	( * distance ) = dist;

	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) y );
	unsigned char * xr;
	xr = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );
	create_rotation ( x, rot, xr );
	
	unsigned char * X;
	unsigned char * Y;

	unsigned int sl = sw . P * ( sw . l ); //section length
	sl = cscmin ( sl, cscmin ( m/2, n/2 ) );

	X = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	Y = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );


	memcpy ( &X[0], &xr[0], sl );
	for ( int i = 0; i < sl; i++ )
		X[sl + i] = DEL;
	memcpy ( &X[sl + sl], &xr[m - sl], sl );
	X[3 * sl] = '\0';
	
	memcpy ( &Y[0], &y[0], sl );
	for ( int i = 0; i < sl; i++ )
		Y[sl + i] = DEL;
	memcpy ( &Y[sl + sl], &y[n - sl], sl );
	Y[3 * sl] = '\0';

	unsigned int mm = sl + sl + sl;
	unsigned int nn = sl + sl + sl;

	double score = -DBL_MAX;
	double max_score = score;
	unsigned int rrot = 0;
	unsigned char * Xr = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	double O = - sw . O;
	double E = - sw . E;
	for ( int i = 0; i < mm; i++ )
	{
		if ( i >= sl && i < 2 * sl )
			continue;
	
		Xr[0] = '\0';
		create_rotation ( X, i, Xr );

		nw ( Xr, mm , Y, nn, O, E, &score, sw . alphabet );	
		if ( score > max_score )
		{
			max_score = score;
			rrot = i;
		}	 
	}
	free ( Xr);

	if ( rrot < sl )
	{
		( * rotation ) = rot + rrot;
	}
	else
	{
		( * rotation ) = rot - ( 3 * sl - rrot );
	}
	( * rotation ) = ( * rotation ) % m;

	free ( xr );
	free ( X );
	free ( Y );
	return EXIT_SUCCESS;
}
