
#ifndef __NCSC_H_INCLUDED__
#define __NCSC_H_INCLUDED__


#include "qcsc.h"
/**
 * hCSC: A heuristic q-gram based approximate pairwise circular sequence alignment algorithm
 */
class nCSC : public qCSC
{
private:
    struct BestMatch runNaive(vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y);

public:
    nCSC(string xx, unsigned int m, string y, unsigned int n, unsigned int q, unsigned int b, string a) : qCSC(xx, m, y, n, q, b, a)
    {
    }
    int run(struct BestMatch * best);
};

#endif
