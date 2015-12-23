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


static struct option long_options[] =
{
   { "method",                  required_argument, NULL, 'm' },
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "q-length",                required_argument, NULL, 'q' },
   { "block-length",            required_argument, NULL, 'l' },
   { "blocks-refine",           optional_argument, NULL, 'P' },
   { "gap-open-penalty",        optional_argument, NULL, 'O' },
   { "gap-extend-penalty",      optional_argument, NULL, 'E' },
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
   sw -> l                              = 10;
   sw -> P                              = 0.0;
   sw -> O                              = 10.0;
   sw -> E                              = 0.5;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "m:a:i:o:q:l:P:O:E:h", long_options, &oi ) ) != - 1 )
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

         case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           args ++;
           break;

         case 'P':
           val = (double) atof ( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           break;

         case 'O':
           val = (double) atof ( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> O = val;
           break;

         case 'E':
           val = (double) atof ( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> E = val;
           break;

         case 'h':
           return ( 0 );
       }
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
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' or `RNA' for nucleotide sequences or\n"
                     "                                      `PROT' for protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename for the rotated sequences.\n" );
   fprintf ( stdout, "  -q, --q-length            <int>     The q-gram length.\n");
   fprintf ( stdout, "  -l, --block-length        <int>     The length of each block.\n");
   fprintf ( stdout, " Extra (Optional and only to be used with saCSC):\n" );
   fprintf ( stdout, "  -P, --blocks-refine       <float>   The number of blocks of length l to use to\n"
                     "                                      refine the results of saCSC by (e.g. 1.0)\n" );
   fprintf ( stdout, "  -O, --gap-open-penalty    <float>   The gap open penalty is the score taken\n"
                     "                                      away when a gap is created. The best\n"
                     "                                      value depends on the choice of comparison\n"
                     "                                      matrix.   The   default   value   assumes\n"
                     "                                      you  are  using  the  EBLOSUM62  matrix\n"
                     "                                      for protein sequences, and the  EDNAFULL\n"
                     "                                      matrix for nucleotide sequences. Floating\n"
                     "                                      point number from 1.0 to 100.0. (default:\n"
                     "                                      10.0)\n" );
   fprintf ( stdout, "  -E, --gap-extend-penalty  <float>   The gap extension penalty is added to\n"
                     "                                      the standard gap penalty for each base or\n"
                     "                                      residue in the gap. This is how long gaps\n"
                     "                                      are penalized. Floating point number from\n"
                     "                                      0.0  to  10.0.  (default:  0.5)\n" );
   fprintf ( stdout, " Other:\n" );
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
