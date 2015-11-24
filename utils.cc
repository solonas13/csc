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

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>
#include "csc.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

//using namespace std;

int EDNA[90];
int BLOSUM[91];

void init_substitution_score_tables ()
{
    int i;
    char edna[] = "ATGCSWRYKMBVHDN";
    for ( i = 0; i < 15; i ++ ) {
	EDNA[(int)edna[i]] = i;
    }
    char blosum[] = "ARNDCQEGHILKMFPSTWYVBZX*";
    for ( i = 0; i < 24; i ++ ) {
	BLOSUM[(int)blosum[i]] = i;
    }
}
 
double delta ( char a, char b, char * alphabet )
{
    /*if ( a == DEL && b == DEL ) {
	return ( double ) 10;
    }
    else if ( ( a == DEL && b != DEL ) || ( b == DEL && a != DEL ) ) {
	return (double) -10;
    }*/

    if ( a == DEL || b == DEL ) {
	a = 'A';
	b = 'A';
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

static struct option long_options[] =
{
   { "method",                  required_argument, NULL, 'm' },
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "q-length-min",            required_argument, NULL, 'q' },
   { "q-length-max",            optional_argument, NULL, 'Q' },
   { "block-length-min",        required_argument, NULL, 'l' },
   { "block-length-max",        optional_argument, NULL, 'L' },
   { "percent-refine",          optional_argument, NULL, 'P' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL, 0   }
};


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
{
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> method                         = NULL;
   sw -> q                              = 5;
   sw -> Q                              = 0;
   sw -> l                              = 10;
   sw -> L                              = 0;
   sw -> P                              = 0.0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "m:a:i:o:q:Q:l:L:P:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
	case 'm':
           sw -> method = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> method, optarg );
           args ++;
           break;
	 
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
           args ++;
           break;

         case 'Q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> Q = val;
           break;

         case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           args ++;
           break;

         case 'L':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> L = val;
           break;

         case 'P':
           val = (double) atof ( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( sw -> Q == 0 )
     {
       sw -> Q = sw -> q;
     }
   if ( sw -> L == 0 )
     {
       sw -> L = sw -> l;
     }

   if ( args < 6 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
}


/* 
Usage of the tool 
*/
void usage ( void )
{
   fprintf ( stdout, " Usage: csc <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -m, --method              <str>     `hCSC' for heuristic, `nCSC' for naive\n"
                     "                                      and `saCSC' for suffix-array algorithm. \n" );
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'\n"
                     "                                      for protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename for the rotated sequences.\n" );
   fprintf ( stdout, "  -q, --q-length-min        <int>     The q-gram length.\n");
   fprintf ( stdout, "  -l, --block-length-min    <int>     The length of each block.\n");
   fprintf ( stdout, " Extra (Optional; Only use with saCSC and hCSC):\n" );
   fprintf ( stdout, "  -Q, --q-length-max        <int>     The maximum q-gram length. The program\n"
                     "                                      will try all in range min .. max.\n" );
   fprintf ( stdout, "  -L, --block-length-max    <int>     The maximum length of each block. The\n"
                     "                                      program will try all in range min .. max.\n" );
   fprintf ( stdout, "  -P, --percent-refine      <float>   Refine the alignment of hCSC/saCSC by\n"
                     "                                      checking a percentage of the ends (e.g. 2.5)\n" );
   fprintf ( stdout, "  -h, --help                <void>    This help message.\n");
}

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
    unsigned int m = strlen ( ( char * ) x );
    memmove ( &rotation[0], &x[offset], m - offset );
    memmove ( &rotation[m - offset], &x[0], offset );
    rotation[m] = '\0';
}

void create_backward_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
    unsigned int m = strlen ( ( char * ) x );
    memmove ( &rotation[0], &x[m - offset], offset );
    memmove ( &rotation[offset], &x[0], m - offset );
    rotation[m] = '\0';
}

int refine ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, double p, char * alphabet )
{
    init_substitution_score_tables();

    //create X and Y prine (Xp, Yp)
    unsigned int sectionLength = ( unsigned int ) floor ( ( p / 100 ) * std::min ( m, n ) );
    unsigned int sl3 = sectionLength * 3;
    unsigned char repeat_char = DEL;
    unsigned char * Xp, * Yp;
    if ( ( Xp = ( unsigned char * ) calloc ( sl3 + 1, sizeof ( unsigned char ) ) ) == NULL ) {
	fprintf ( stderr, " Error: could not allocate Xp.\n" );
	return 1;
    }
    if ( ( Yp = ( unsigned char * ) calloc ( sl3 + 1, sizeof ( unsigned char ) ) ) == NULL ) {
	fprintf ( stderr, " Error: could not allocate Yp.\n" );
	return 1;
    }
    
    //make Xp and Yp contain prefix, repeat_char middle and suffix of x and y respectively, e.g. AATGCA$$$$$GGGAT
    memcpy ( Xp, x, sectionLength * sizeof ( unsigned char ) );
    memset ( Xp + sectionLength * sizeof ( unsigned char ), repeat_char, sectionLength * sizeof ( unsigned char ) );
    memcpy ( Xp + 2 * sectionLength * sizeof ( unsigned char ), &x[ m - sectionLength ], sectionLength * sizeof ( unsigned char ) );
    Xp[3 * sectionLength] = '\0';
    memcpy ( Yp, y, sectionLength * sizeof ( unsigned char ) );
    memset ( Yp + sectionLength * sizeof ( unsigned char ), repeat_char, sectionLength * sizeof ( unsigned char ) );
    memcpy ( Yp + 2 * sectionLength * sizeof ( unsigned char ), &y[ n - sectionLength ], sectionLength * sizeof ( unsigned char ) );
    Yp[3 * sectionLength] = '\0';
    
    unsigned int i, j, r, rotation;
    double max_score = -DBL_MAX;
    double g = -5; //open gap penalty
    double h = -0.5; //extend gap penalty

    double * d0;
    double * d1;
    double * t0;
    double * t1;
    double * in;
    if ( ( d0 = ( double * ) calloc ( sl3 + 1 , sizeof ( double ) ) ) == NULL )
    {
	fprintf ( stderr, " Error: 'd0' could not be allocated!\n");
	return 0;
    }
    if ( ( d1 = ( double * ) calloc ( sl3 + 1 , sizeof ( double ) ) ) == NULL  )
    {
	fprintf ( stderr, " Error: 'd1' could not be allocated!\n");
	return 0;
    }
    if ( ( t0 = ( double * ) calloc ( sl3 + 1 , sizeof ( double ) ) ) == NULL )
    {
	fprintf ( stderr, " Error: 't0' could not be allocated!\n");
	return 0;
    }
    if ( ( t1 = ( double * ) calloc ( sl3 + 1 , sizeof ( double ) ) ) == NULL )
    {
	fprintf ( stderr, " Error: 't1' could not be allocated!\n");
	return 0;
    }
    if ( ( in = ( double * ) calloc ( sl3 + 1 , sizeof ( double ) ) ) == NULL )
    {
	fprintf ( stderr, " Error: 'in' could not be allocated!\n");
	return 0;
    }

    unsigned char * yr;
    if ( ( yr = ( unsigned char * ) calloc ( ( sl3 + 1 ) , sizeof ( unsigned char ) ) ) == NULL )
    {
	fprintf( stderr, " Error: 'yr' could not be allocated!\n");
	return 0;
    }

    for ( r = 0; r < sl3; r++ )
    {
        if ( r >= sectionLength && r < 2 * sectionLength ) {
	    continue;
	}

	yr[0] = '\0';

	create_rotation ( Xp, r, yr );

	memset ( d0, 0, sizeof ( double ) * sl3 + 1 );
	memset ( d1, 0, sizeof ( double ) * sl3 + 1 );
	memset ( t0, 0, sizeof ( double ) * sl3 + 1 );
	memset ( t1, 0, sizeof ( double ) * sl3 + 1 );
	memset ( in, 0, sizeof ( double ) * sl3 + 1 );

	for ( j = 0; j < sl3 + 1; j++ )
	{
	    d0[j] = -DBL_MAX;
	    in[j] = -DBL_MAX;
	}

	t0[0] = 0;
	t0[1] = g;
	t1[0] = g; 

	for ( j = 2; j < sl3 + 1; j++ ) {
	    t0[j] = t0[j - 1] + h;
	}

	for ( i = 1; i < sl3 + 1; i++ )
	{
	    for ( j = 0; j < sl3 + 1; j++ )
	    {
		double u, v, w;

		switch ( i % 2 ) 
		{

		  case 0:

		    if ( j == 0 )
		    {
			d0[j] = -DBL_MAX;
			in[j] = -DBL_MAX;
			if ( i >= 2 ) {
			    t0[0] = t1[0] + h;
			}
		    }
		    else 
		    {
			d0[j] = std::max ( d1[j] + h, t1[j] + g );
			u = d0[j];

			in[j] = std::max ( in[j - 1] + h, t0[j - 1] + g ); //i0
			v = in[j];

			w = t1[j - 1] + delta ( Yp[j - 1], yr[i - 1], alphabet );

			t0[j] = std::max ( w, std::max ( u, v ) );

			if ( i == sl3 && j == sl3 && t0[j] > max_score )
			{
			    max_score = t0[j];
			    rotation  = ( r >= sectionLength ) ? -( sl3 - r ) : r;
			}
		    }

		    break;

		  case 1:

		    if ( j == 0 )
		    {
			d1[j] = -DBL_MAX;
			in[j] = -DBL_MAX;
			if ( i >= 2 ) {
			    t1[0] = t0[0] + h;
			}
		    }	
		    else 
		    {
			d1[j] = std::max ( d0[j] + h, t0[j] + g );
			u = d1[j];

			in[j] = std::max ( in[j - 1] + h, t1[j - 1] + g ); //i1
			v = in[j];

			w = t0[j - 1] + delta (  Yp[j - 1], yr[i - 1], alphabet );

			t1[j] = std::max ( w, std::max ( u, v ) );

			if ( i == sl3 && j == sl3 && t1[j] > max_score )
			{
			    max_score = t1[j];
			    rotation  = ( r >= sectionLength ) ? -( sl3 - r ) : r;
			}
		    }

		    break;

		}

	    }

	}

    }

    free ( Xp );
    free ( Yp );
    free ( yr );
    free( d0 );
    free( d1 );
    free( t0 );
    free( t1 );
    free( in );

    return rotation;
}
