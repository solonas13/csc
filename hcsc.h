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

#ifndef __HCSC_H_INCLUDED__
#define __HCSC_H_INCLUDED__

#include "qcsc.h"

/**
 * hCSC: A heuristic q-gram based approximate pairwise circular sequence alignment algorithm
 */
class hCSC : public qCSC
{
private:
    struct BestMatch getBestScoringBlock(vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y);
    struct BestMatch refine(struct BestMatch oldBest, vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y);

public:
    hCSC(string xx, unsigned int m, string y, unsigned int n, unsigned int q, unsigned int b, string a) : qCSC(xx, m, y, n, q, b, a){}
    int run(struct BestMatch * best);
};

#endif
