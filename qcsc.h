
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

#ifndef __QCSC_H_INCLUDED__
#define __QCSC_H_INCLUDED__


#include <cmath>
#include <climits>
#include <cstring>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

typedef unsigned long int       WORD;
#define WORD_SIZE               sizeof(WORD) * 8

/**
 * A struct to hold the best matching block information -
 * the block index, score and block starting position for
 * the naive and heuristic methods
 */
struct BestMatch {
    unsigned int index;
    unsigned int score;
    unsigned int pos;
};

/**
 * qCSC: An abstract class that is extended by hCSC and nCSC to do their q-gram based pairwise circular sequence alignment
 * 
 * @abstract
 */
class qCSC
{
protected:

    /**
     * @var sigma Alphabet size
     */
    unsigned int sigma;
    /**
     * @var xx The long text string (doubled-up x)
     */
    string xx;
    /**
     * @var m Length of xx
     */
    unsigned int m;
    /**
     * @var y The search string to find in xx
     */
    string y;
    /**
     * @var n Length of y
     */
    unsigned int n;
    /**
     * @var bxSize The size of a block of x in characters
     */
    unsigned int bxSize;
    /**
     * @var bxNum The number of blocks in xx
     */
    unsigned int bxNum;
    /**
     * @var qxNum The number of qgrams in an x block
     */
    unsigned int qxNum;
    /**
     * @var bySize The size of a block of y in characters
     */
    unsigned int bySize;
    /**
     * @var qyNum The number of qgrams in a y block
     */
    unsigned int qyNum;
    /**
     * @var byNum The number of blocks in y
     */
    unsigned int byNum;
    /**
     * @var l How many bits are required to represent a character of the Alphabet
     */
    unsigned int l;
    /**
     * @var qSize The size of the q-gram in characters
     */
    unsigned int qSize;
    /**
     * @var chars Stores the decimal value of characters of the alphabet
     */
    WORD * chars;
    /**
     * @var numPerms The number of permutations possible (sigma^q)
     */
    unsigned int numPerms;

    /**
     * Fills chars array and returns the size of the alphabet
     *
     * @return Number of unique characters
     */
    unsigned int generateCharHashmap(void)
    {
	char * c;
	if ((c = (char *) calloc(CHAR_MAX, sizeof(char))) == NULL) {
	    cerr << "Error assigning char" << endl;
	}
      
	//go through xx and y and mark characters found
	unsigned int i;
	for (i = 0; i < this->m / 2; i++) {
	    c[(char)this->xx[i]] = 1;
	}
	for (i = 0; i < this->n; i++) {
	    c[(char)this->y[i]] = 1;
	}

	//give characters unique id
	WORD counter = 0;
	for (i = 0; i < CHAR_MAX; i++) {
	    if (c[i] == 1) {
		//cout << (char)i << counter << endl;
		this->chars[i] = counter++;
	    }
	}

	free(c);

	return (unsigned int)counter;
    }

    /**
     * Turns string's characters into Q-Grams and shuffles them onto a word
     * Warning: No word-size bounds checking
     *
     * @param s A string
     * @return
     */
    WORD shuffleOntoWord(string s)
    {
	WORD w = 0;

	unsigned int i;
	for (i = 0; i < s.length(); i++) {
	    w = w << this->l;
	    w = w | this->chars[s.at(i)];
	}

	return w;
    }

    /**
     * Fills a matrix with the calculated qgrams
     * 
     * @param Z The matrix to fill
     * @param s The sequence to put into the matrix
     * @param len The length of s
     * @param charsInBlock The number of characters in a block of s
     * @param qgramsInBlock The number of qGrams in a block of s
     * @return
     */
    vector<unordered_map<WORD, unsigned int>> fillQGramBlocks(vector<unordered_map<WORD, unsigned int>> Z, string s, unsigned int len, unsigned int charsInBlock, unsigned int qgramsInBlock)
    {
	int i, j, blkLen, blockNum = 0, qSize = (int)this->qSize;
	WORD qGram;

	//loop through string in block-sized segments
	for (i = 0; i < len; i += charsInBlock) {

	    //get the current block
	    blkLen = s.substr(i, qgramsInBlock).length();

	    //count q-grams in each block
	    for (j = 0; j < blkLen - qSize + 1; j++) {
		//read the current qGram from the block
		qGram = shuffleOntoWord(s.substr(i + j, qSize));

		//initialise q-gram count or increment if already exists
		Z[blockNum].emplace(qGram, 0);
		Z[blockNum][qGram]++;
	    }

	    //next block index
	    blockNum++;
	}

	return Z;
    }

    /**
     * Finds a value in a map and returns it or returns 0. This method does not
     * create a new element/return pointer as default C++ behaviour does so better
     * 
     * @param map
     * @param key
     * @return value stored for the given key or 0
     */
    int findInMap(unordered_map<WORD, unsigned int> map, WORD key) {
	unordered_map<WORD, unsigned int>::const_iterator i = map.find(key);
	return (i == map.end()) ? 0 : i->second;
    }

    /**
     * Finds a value in a map and returns it or returns 0. This method does not
     * create a new element/return pointer as default C++ behaviour does so better
     * 
     * @param map
     * @param key
     * @return value stored for the given key or 0
     */
    int findInMap(unordered_map<WORD, int> map, WORD key) {
	unordered_map<WORD, int>::const_iterator i = map.find(key);
	return (i == map.end()) ? 0 : i->second;
    }

public:

    /**
     * The hCSC class constructor
     *
     * @param xx long text string (you need to pass it in doubled if required)
     * @param m length of xx
     * @param y search string to find in xx
     * @param n length of y
     * @param q q-gram length as set by the user
     * @param b block length as set by the user
     * @param a alphabet characters as defined by the user
     * @return
     */
    qCSC(string xx, unsigned int m, string y, unsigned int n, unsigned int q, unsigned int b, string a)
    {
	this->xx = xx;
	this->m = m;
	this->y = y;
	this->n = n;
	
	//initialise the chars array - faster than hashmap look-up but takes up a little bit more memory
	if ((this->chars = (WORD *) calloc(CHAR_MAX, sizeof(WORD))) == NULL) {
	    cerr << "Could not allocate the chars map." << endl;
	}

	//sigma holds the size of the alphabet, chars the characters
	if (a.length() == 0) {
	    this->sigma = this->generateCharHashmap();
	} else {
	    this->sigma = a.length();
	    WORD i;
	    for (i = 0; i < this->sigma; i++) {
		this->chars[a[i]] = i;
	    }
	}
	//cout << "sigma: " << this->sigma << endl;

	//l holds the size of a single letter in the q-gram in bits; If sigma is 4, l is 2 (00,01,10,11)
	this->l = (unsigned int) ceil(sqrt(this->sigma));
	//cout << "l: " << this->l << endl;

	//qSize holds the q-gram size
	if (q == 0) {
	    //qx holds the q-gram size of x; If x is 250, sigma is 4, q = 4
	    unsigned int qx = (unsigned int) ceil((double)log(0.5 * this->m) / (double)log(this->sigma));
	    //qy holds the q-gram size of y; If y is 250, sigma is 4, q = 4
	    unsigned int qy = (unsigned int) ceil((double)log(this->n) / (double)log(this->sigma));
	    //qgram set to minimum q-gram size of the two
	    this->qSize = min(qx, qy);
	} else {
	    //set to user defined size
	    this->qSize = q;
	}
	//cout << "Q-Gram size: " << this->qSize << endl;

	//numPerms is the number of q-sized permutations that can be made from the alphabet
	this->numPerms = (unsigned int) pow((double)this->sigma, (double)this->qSize);
	//cout << "perms: " << this->numPerms << endl;

	//how many characters in an xx block
	if (b == 0) {
	    this->bxSize = (unsigned int) ceil((double)sqrt(0.5 * this->m));
	} else {
	    this->bxSize = b;
	}
	//cout << "XX Block size: " << this->bxSize << endl;

	//how many qgrams in an x block
	this->qxNum = this->bxSize + this->qSize - 1;
	//cout << "num q in x: " << this->qxNum << endl;

	//how many blocks in xx
	this->bxNum = (unsigned int) ceil((double)this->m / (double)this->bxSize);
	//cout << "blocks x: " << this->bxNum << endl;

	//how many characters in an y block
	if (b == 0) {
	    this->bySize = (unsigned int) ceil((double)sqrt(this->n));
	} else {
	    this->bySize = b;
	}
	//cout << "Y Block size: " << this->bySize << endl;

	//how many qgrams in an y block
	this->qyNum = this->bySize + this->qSize - 1;
	//cout << "num q in yb: " << this->qyNum << endl;

	//how many blocks in y
	this->byNum = (unsigned int) ceil((double)this->n / (double)this->bySize);
	//cout << "blocks y: " << this->byNum << endl;
    }

    /**
     * Destructor to free reserved memory
     */
    ~qCSC()
    {
	free(this->chars);
    }

    /**
     * @abstract
     */
    virtual int run(struct BestMatch * best) = 0;
};

#endif

