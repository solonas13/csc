
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

#ifndef __CYC_NW__
#define __CYC_NW__

#include <iostream>
#include <algorithm>
#include <cfloat>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>
#include "EDNAFULL.h"

#define ALLOC_SIZE              1048576
#define ALPHABET                "ATGCN"

struct option long_options[] =
{
  { "input-file",              required_argument, NULL, 'i' },
  { "output-file",             required_argument, NULL, 'o' },
  { "open-gap-penalty",        required_argument, NULL, 'O' },
  { "extend-gap-penalty",      required_argument, NULL, 'E' },
  { "help",                    no_argument,       NULL, 'h' },
  { NULL,                      0,                 NULL, 0   }
};

struct TSwitch {
    char * input_filename;
    char * output_filename;
    double open_gap_penalty;
    double extend_gap_penalty;
};

using namespace std;

unsigned int nuc_char_to_index ( char a );
int nuc_delta ( char a, char b );
unsigned int create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
unsigned int cyc_nw_ls ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, double o, double e, double * score, int * rot );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
double gettime( void );

#endif

