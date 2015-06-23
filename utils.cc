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
#include <sys/time.h>
#include "csc.h"

static struct option long_options[] =
 {
   { "method",                  required_argument, NULL, 'm' },
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "q-length-min",            required_argument, NULL, 'q' },
   { "q-length-max",            required_argument, NULL, 'Q' },
   { "block-length-min",        required_argument, NULL, 'l' },
   { "block-length-max",        required_argument, NULL, 'L' },
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
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "m:a:i:o:q:Q:l:L:h", long_options, &oi ) ) != - 1 )
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

         case 'L':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> L = val;
           args ++;
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
   fprintf ( stdout, "  -l, --block-length-max    <int>     The maximum length of each block. The\n"
                     "                                      program will try all in range min .. max.\n" );
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
