/**
    CSC: Circular Sequence Comparison
    Copyright (C) 2015 Solon P. Pissis. 

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
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "csc.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                  // the input file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char ** seq    = NULL;         	// the sequence in memory
        unsigned char ** seq_id = NULL;         	// the sequence id in memory
	char *          alphabet;               // the alphabet
	unsigned int    i, j;
	unsigned int    b, q;	

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 9 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet = ( char * ) DNA;
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

		q       = sw . q;
		b       = sw . b;

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

	if ( sw . b < 1 || sw . b > m - sw . q + 1  || sw . b > n - sw . q + 1 )
	{
        	fprintf( stderr, " Error: Illegal number of blocks.\n" );
		return ( 1 );
	}

	if ( sw . q >= m || sw . q >= n )
	{
        	fprintf( stderr, " Error: Illegal q-gram length.\n" );
		return ( 1 );
	}

	if ( num_seqs > 2 )
	{
        	fprintf( stderr, " Warning: %d sequences were read from file %s.\n", num_seqs, input_filename );
        	fprintf( stderr, " Warning: Only the first two (%s, %s) will be processed!\n", seq_id[0], seq_id[1] );
	}
	
	unsigned int distance = m + n;
	unsigned int rotation = 0;

	#if 0
	for ( sw . b = 500; sw . b <= 1500; sw . b++ )
	{	
		circular_sequence_comparison ( seq[0], seq[1], sw, &rotation, &distance );
	}
	#endif

	circular_sequence_comparison ( seq[0], seq[1], sw, &rotation, &distance );

	TPOcc D;
	D . err = distance;
	D . rot = rotation;

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

	double end = gettime();

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	unsigned char * rot_str = ( unsigned char * ) calloc ( m + 1, sizeof ( unsigned char ) );
	create_rotation ( seq[0], D . rot, rot_str );
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
        fprintf( stderr, " q-gram length is %d\n", sw . q );
        fprintf( stderr, " Number of blocks is %d\n", sw . b );
        fprintf( stderr, " Block length is %d\n", m / sw . b );
        fprintf( stderr, " Blockwise q-gram distance: %d\n", ( int ) D . err );
        fprintf( stderr, " Rotation                 : %d\n", D . rot );
        fprintf( stderr, " (Multi)FASTA output file : %s\n", sw . output_filename );
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
