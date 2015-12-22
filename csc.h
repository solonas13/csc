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

#ifndef __CSC__
#define __CSC__

#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"
#define METHOD_H                "hCSC"
#define METHOD_N                "nCSC"
#define METHOD_SA               "saCSC"
#define ALPHABET_DNA            "DNA"
#define ALPHABET_RNA            "RNA"
#define ALPHABET_PROT           "PROT"
#define ALPHABET_IUPAC          "IUPAC"
#define DNA                     "ACGTN"                         //DNA alphabet
#define RNA                     "ACGUN"                         //RNA alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC nucleotide alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYVX"         //Proteins alphabet

#define cscmax(a,b) ((a) > (b)) ? (a) : (b)
#define cscmin(a,b) ((a) < (b)) ? (a) : (b)

struct TSwitch
{
    char *               input_filename;         // the input file name
    char *               output_filename;        // the output file name
    char *               alphabet;               // the output file name
    char *               method;                 // algorithm/method
    unsigned int         l;                      // block length
    unsigned int         q;                      // q-gram size
    double               P;                      // (optional) number of blocks to refine by
    double               O;
    double               E;
};

struct TPOcc
{
    unsigned int         err;
    unsigned int         rot;
};

extern int EDNA[];
extern int BLOSUM[];

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
void create_backward_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
int refine ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, double p, char * alphabet );
void init_substitution_score_tables ();
double delta ( char a, char b, char * alphabet );

#endif
