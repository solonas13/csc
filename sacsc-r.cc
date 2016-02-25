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

unsigned int nw ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, double o, double e, double * score, char * alphabet, double * d0, double * d1, double * t0, double * t1, double * in )
{
	int i, j;
	double g = o;
	double h = e;
	double u, v, w;
	double max_score = -DBL_MAX;
	
        for ( j = 0; j < n + 1; j++ )
	{
		d1[j] = 0;
		t1[j] = 0;
		d0[j] = -DBL_MAX;
		in[j] = -DBL_MAX;
	}

	if ( m > 0 )	t1[0] = g; 
	if ( n > 0 )	t0[1] = g;
	t0[0] = 0;

	for ( j = 2; j < n + 1; j++ )
		t0[j] = t0[j - 1] + h;

	for( i = 1; i < m + 1; i++ )
	{
        	for( j = 0; j < n + 1; j++ )
        	{
			switch ( i % 2 ) 
			{
			    case 0:

				if ( j == 0 )
				{
					d0[j] = -DBL_MAX;
					in[j] = -DBL_MAX;
					if ( i >= 2 )
						t0[0] = t1[0] + h;
				}
				else 
				{
					u = d0[j] = cscmax ( d1[j] + h, t1[j] + g );

					v = in[j] = cscmax ( in[j - 1] + h, t0[j - 1] + g ); //i0

					w = t1[j - 1] + delta( p[j - 1], t[i - 1], alphabet );

					t0[j] = cscmax ( w, cscmax ( u, v ) );

					if ( i == m && j == n && t0[n] > max_score )
					{
						max_score = t0[j];
					}
				}

				break;
				
			    case 1:

				if ( j == 0 )
				{
					d1[j] = -DBL_MAX;
					in[j] = -DBL_MAX;
					if ( i >= 2 )
						t1[0] = t0[0] + h;
				}	
				else 
				{
					u = d1[j] = cscmax ( d0[j] + h, t0[j] + g );

					v = in[j] = cscmax ( in[j - 1] + h, t1[j - 1] + g ); //i1

					w = t0[j - 1] + delta( p[j - 1], t[i - 1], alphabet );

					t1[j] = cscmax ( w, cscmax ( u, v ) );

					if ( i == m && j == n && t1[n] > max_score )
					{
						max_score = t1[j];
					}
				}

				break;
			}
        	}
    	}

	( * score ) = max_score;

	return EXIT_SUCCESS;
}

unsigned int sacsc_refinement (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance )
{
	unsigned int i;
	unsigned int rot;
	unsigned int dist;
	double starttime = gettime();
	
	circular_sequence_comparison ( x, y, sw, &rot, &dist );

	double endtime = gettime() - starttime;

	fprintf ( stderr, "saCSC completed in %fs with initial rotation %u.\n", endtime, rot );

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

	init_substitution_score_tables ();
	
	double * d0;
	double * d1;
	double * t0;
	double * t1;
	double * in;
	if ( ( d0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
	{
	    fprintf( stderr, " Error: 'd0' could not be allocated!\n");
	    return EXIT_FAILURE;
	}
	if ( ( d1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL  )
	{
	    fprintf( stderr, " Error: 'd1' could not be allocated!\n");
	    return EXIT_FAILURE;
	}
	if ( ( t0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
	{
	    fprintf( stderr, " Error: 't0' could not be allocated!\n");
	    return EXIT_FAILURE;
	}
	if ( ( t1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
	{
	    fprintf( stderr, " Error: 't1' could not be allocated!\n");
	    return EXIT_FAILURE;
	}
	if ( ( in = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
	{
	    fprintf( stderr, " Error: 'in' could not be allocated!\n");
	    return EXIT_FAILURE;
	}

	fprintf ( stderr, "Refining using a sequence length of %u...\n", sl );

	for ( int i = 0; i < mm; i++ )
	{
		if ( i >= sl && i < 2 * sl )
			continue;
	
		Xr[0] = '\0';
		create_rotation ( X, i, Xr );

		nw ( Xr, mm , Y, nn, O, E, &score, sw . alphabet, d0, d1, t0, t1, in );	
		if ( score > max_score )
		{
			max_score = score;
			rrot = i;
		}	 
	}
	free ( Xr );

	free( d0 );
	free( d1 );
	free( t0 );
	free( t1 );
	free( in );

	int final_rot;
	if ( rrot < sl )
	{
		final_rot = rot + rrot;
	}
	else
	{
		final_rot = rot - ( 3 * sl - rrot );
	}

	if ( final_rot > ( int ) m )
	{
		( * rotation ) = final_rot % m;	
	}
	else if ( final_rot < 0 )
	{
		( * rotation ) = m + final_rot;
	}
	else
		( * rotation ) = final_rot;

	free ( xr );
	free ( X );
	free ( Y );
	return EXIT_SUCCESS;
}
