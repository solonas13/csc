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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/time.h>
#include "csc.h"
#include "sacsc.h"
#include "hcsc.h"
#include "ncsc.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *           in_fd;                  // the input file descriptor
	FILE *           out_fd;                 // the input file descriptor
        char *           input_filename;         // the input file name
        char *           output_filename;        // the output file name
        char *           method;                 // the method - nCSC, hCSC, saCSC
        unsigned char ** seq    = NULL;          // the sequence in memory
        unsigned char ** seq_id = NULL;          // the sequence id in memory
	char *           alphabet;               // the alphabet
	unsigned int     l, L, q, Q;             // the program parameters
	double           P;                      // the program parameters
	
	unsigned int     h, i, j, k;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 10 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( METHOD_H, sw . method ) )   method = ( char * ) METHOD_H;
                else if ( ! strcmp ( METHOD_N, sw . method ) )   method = ( char * ) METHOD_N;
                else if ( ! strcmp ( METHOD_SA, sw . method ) )  method = ( char * ) METHOD_SA;
                else
                {
                        fprintf ( stderr, " Error: Method argument should be `hCSC', `nCSC' or `saCSC' for heuristic, naive or suffix-array Circular Sequence Comparison.\n" );
                        return ( 1 );
                }

                if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet = ( char * ) DNA;
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

                if      ( sw . L < sw . l || sw . Q < sw . q )
                {
                        fprintf ( stderr, " Error: The optional maxmimum parameters cannot be less than the mandatory parameters!\n" );
                        return ( 1 );
                }

                if      ( sw . P < 0 || sw . P >= 50.0 )
                {
                        fprintf ( stderr, " Error: The optional percentage flag should be in the range of 0 to 50%%.\n" );
                        return ( 1 );
                }

                q       = sw . q;
                Q       = sw . Q;
                l       = sw . l;
                L       = sw . L;
		P       = sw . P;

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;
        }

	double start = gettime();

        /* Read the (Multi)FASTA file in memory */
        fprintf ( stderr, " Reading the (Multi)FASTA input file: %s\n", input_filename );
        if ( ! ( in_fd = fopen ( input_filename, "r") ) )
        {
                fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
                return ( 1 );
        }

        char c;
        unsigned int num_seqs = 0;           	// the total number of sequences considered
        unsigned int total_length = 0;          // the total number of sequences considered
        unsigned int max_alloc_seq_id = 0;
        unsigned int max_alloc_seq = 0;
        c = fgetc( in_fd );
        do
        {
                if ( c != '>' )
                {
                        fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
                        return ( 1 );
                }
                else
                {
                        if ( num_seqs >= max_alloc_seq_id )
                        {
                                seq_id = ( unsigned char ** ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
                                max_alloc_seq_id += ALLOC_SIZE;
                        }

                        unsigned int max_alloc_seq_id_len = 0;
                        unsigned int seq_id_len = 0;

                        seq_id[ num_seqs ] = NULL;

                        while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
                        {
                                if ( seq_id_len >= max_alloc_seq_id_len )
                                {
                                        seq_id[ num_seqs ] = ( unsigned char * ) realloc ( seq_id[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                        max_alloc_seq_id_len += ALLOC_SIZE;
                                }
                                seq_id[ num_seqs ][ seq_id_len++ ] = c;
                        }
                        seq_id[ num_seqs ][ seq_id_len ] = '\0';

                }
		if ( num_seqs >= max_alloc_seq )
                {
                        seq = ( unsigned char ** ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
                        max_alloc_seq += ALLOC_SIZE;
                }

                unsigned int seq_len = 0;
                unsigned int max_alloc_seq_len = 0;

                seq[ num_seqs ] = NULL;

                while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
                {
                        if( seq_len == 0 && c == '\n' )
                        {
                                fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
                                c = fgetc( in_fd );
                                break;
                        }
                        if( c == '\n' || c == ' ' ) continue;

                        c = toupper( c );

                        if ( seq_len >= max_alloc_seq_len )
                        {
                                seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                max_alloc_seq_len += ALLOC_SIZE;
                        }

                        if( strchr ( alphabet, c ) )
                        {
                                seq[ num_seqs ][ seq_len++ ] = c;
                        }
                        else
                        {
                                fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
                                return ( 1 );
                        }

                }

                if( seq_len != 0 )
                {
                        if ( seq_len >= max_alloc_seq_len )
                        {
                                seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                max_alloc_seq_len += ALLOC_SIZE;
                        }
                        seq[ num_seqs ][ seq_len ] = '\0';
                        total_length += seq_len;
                        num_seqs++;
                }

        } while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	unsigned int m = strlen ( ( char * ) seq[0] );
	unsigned int n = strlen ( ( char * ) seq[1] );

	if ( sw . L < 1 || sw . L > m - sw . q + 1  || sw . L > n - sw . q + 1 )
	{
        	fprintf( stderr, " Error: Illegal block length.\n" );
		return ( 1 );
	}

	if ( sw . Q >= sw . l )
	{
        	fprintf( stderr, " Error: Illegal q-gram length.\n" );
		return ( 1 );
	}

	if ( num_seqs > 2 )
	{
        	fprintf( stderr, " Warning: %d sequences were read from file %s.\n", num_seqs, input_filename );
        	fprintf( stderr, " Warning: Only the first two (%s, %s) will be processed!\n", seq_id[0], seq_id[1] );
	}


	/* Run the algorithm using the user's chosen method - Naive method does not go through repeats */
	string xx ( ( char * ) seq[0] );
	xx = xx + xx;
	string y ( ( char * ) seq[1] );
	string alphabetLetters = "";
	if ( strcmp ( sw . alphabet, ALPHABET_DNA ) == 0 )
	{
		alphabetLetters = DNA;
	}
	else if ( strcmp ( sw. alphabet, ALPHABET_PROT ) == 0 )
	{
		alphabetLetters = PROT;
	}
	else
	{
		alphabetLetters = IUPAC;
	}
	unsigned int distance = m + n;
	unsigned int rotation = 0;
	TPOcc D, DTemp;
	D . err = UINT_MAX;
	struct BestMatch bm;
	struct TSwitch swTemp;
	memcpy ( &swTemp, &sw, sizeof ( TSwitch ) );

	if ( strcmp ( method, METHOD_N ) == 0 )
	{
		nCSC ncsc ( xx, 2 * m, y, n, sw . q, sw . l, alphabetLetters );
		ncsc . run ( &bm );
		D . err = bm . score;
		D . rot = bm . pos;
	}
	else
	{
		unsigned int q, l;
		for ( q = sw . q ; q <= sw . Q; q ++ ) {
			for ( l = sw . l ; l <= sw . L; l ++ ) {
				swTemp . q = q;
				swTemp . l = l;
			  
				if ( strcmp ( method, METHOD_SA ) == 0 )
				{
				        swTemp . l = (unsigned int) ( m / l ); //use block number instead of length for saCSC				  
					circular_sequence_comparison ( seq[0], seq[1], swTemp, &rotation, &distance );
					DTemp . err = distance;
					DTemp . rot = rotation;
				} else {
					hCSC hcsc ( xx, 2 * m, y, n, sw . q, sw . l, alphabetLetters );
					hcsc . run ( &bm );
					DTemp . err = bm . score;
					DTemp . rot = bm . pos;
				}

				//store best combination
				if ( DTemp . err < D . err ) {
				    D . err = DTemp . err;
				    D . rot = DTemp . rot;
				    sw . q = q;
				    sw . l = l;
				}
			}
		}
	}

	
	#if 0
	for ( int i = 0; i < num_seqs; i++ )
	{
		unsigned int m = strlen ( ( char * ) seq[i] );
		for ( int j = 0; j < num_seqs; j ++ )
		{
			if ( i == j ) continue;

			unsigned int n = strlen ( ( char * ) seq[j] );

			/* Initialise the arrays */
			D[i][j] . err = n + m;
			unsigned int distance = n + m;
			unsigned int rotation = 0;

			circular_sequence_comparison ( seq[i], seq[j], sw, &rotation, &distance );

			D[i][j] . err = distance;
			D[i][j] . rot = rotation;
		}
	}
	#endif

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	unsigned char * rot_str;
	if ( ( rot_str = ( unsigned char * ) calloc ( m + 1, sizeof ( unsigned char ) ) ) == NULL )
	{
		fprintf ( stderr, " Error: Could not allocate rot_str!\n" );
		return ( 1 );
	}

	if ( strcmp ( method, METHOD_N ) == 0 || sw . P == 0 )
	{
		create_rotation ( seq[0], D . rot, rot_str );
	}
	//refine the result if P is defined for method hCSC or saCSC
	else if ( sw . P > 0 && strcmp ( method, METHOD_H ) == 0 || strcmp ( method, METHOD_SA ) == 0 )
	{
		create_rotation ( seq[0], D . rot, rot_str );
		memcpy ( seq[0], rot_str, m * sizeof ( unsigned char ) );
		rot_str[0] = '\0';
		int refinement = refine ( seq[0], m, seq[1], n, sw . P, sw . alphabet );

		if ( refinement >= 0 ) {
		    create_rotation ( seq[0], refinement, rot_str );
		} else {
		    create_backward_rotation ( seq[0], -refinement, rot_str );
		}

		D . rot = (unsigned int) (int)D . rot + refinement;
	}


	double end = gettime();


	fprintf( out_fd, ">%s\n", seq_id[0] );
	fprintf( out_fd, "%s\n", rot_str );
	free ( rot_str );
	fprintf( out_fd, ">%s\n", seq_id[1] );
	fprintf( out_fd, "%s\n", seq[1] );


	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

        fprintf( stderr, " Seq x id is %s and its length is %d\n", seq_id[0], m );
        fprintf( stderr, " Seq y id is %s and its length is %d\n", seq_id[1], n );
        fprintf( stderr, " q-gram length is %d\n",                 sw . q );
        fprintf( stderr, " Number of blocks is %d\n",              m / sw . l );
        fprintf( stderr, " Block length is %d\n",                  sw . l );
        fprintf( stderr, " Blockwise q-gram distance: %u\n",       D . err );
        fprintf( stderr, " Rotation                 : %u\n",       D . rot );
        fprintf( stderr, " (Multi)FASTA output file : %s\n",       sw . output_filename );
        fprintf( stderr, "Elapsed time for comparing sequences: %lf secs\n", ( end - start ) );

	/* De-allocate */
        for ( i = 0; i < num_seqs; i ++ )
        {
                free ( seq[i] );
                free ( seq_id[i] );
        }
        free ( seq );
        free ( seq_id );

        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );

	return ( 0 );
}
