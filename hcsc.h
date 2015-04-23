
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
    hCSC(string xx, unsigned int m, string y, unsigned int n, unsigned int q, unsigned int b, string a) : qCSC(xx, m, y, n, q, b, a)
    {
    }
    
    int run(struct BestMatch * best);
};

#endif
