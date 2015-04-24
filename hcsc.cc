
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

#include "hcsc.h"

/////////////////////
// Private Methods //
/////////////////////

/**
  * Compares the xx string q-gram signature against y's q-grams and returns
  * the starting block of xx with the best score
  * 
  * @param XX q-grams block matrix
  * @param Y q-grams block matrix
  * @return
  */
struct BestMatch hCSC::getBestScoringBlock(vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y)
{
    unsigned int numBlocks = 1 + this->bxNum - this->byNum;
    unsigned int score, tempScore, bestScoreIndex, lowestScore = UINT_MAX;
    unsigned int i = 0, j, k;

    unordered_map<WORD, int> x;
    x.reserve(this->qyNum + this->qxNum);

    //loop through blocks
    do {

	//find best score
	score = 0;
	for (j = i, k = 0; k < this->byNum; j++, k++) {
	    tempScore = 0;

	    //get qgrams scores for y block
	    for (auto py = Y[k].begin(); py != Y[k].end(); py++) {
		tempScore += abs(this->findInMap(XX[j], py->first) - this->findInMap(Y[k], py->first));
		x[py->first] = -1;
	    }

	    //get qgrams scores for x block
	    for (auto px = XX[j].begin(); px != XX[j].end(); px++) {
		if (this->findInMap(x, px->first) != -1) {
		    tempScore += abs(this->findInMap(XX[j], px->first) - this->findInMap(Y[k], px->first));
		}
	    }

	    score += tempScore;

	    x.clear();
	}

	//record best score
	if (score < lowestScore) {
	    lowestScore = score;
	    bestScoreIndex = i;
	}

	i++;
    }
    while (i < numBlocks);

    struct BestMatch best;
    best.index = bestScoreIndex; //block number starting from 0
    best.score = lowestScore;
    best.pos = this->bxSize * bestScoreIndex; //index on xx where the block starts

    return best;
}

/**
  * Refine takes the best match found using the heuristic method as a
  * parameter and looks m/k-1 places either side of the suggested block.
  * It removes the q-gram left of the previous block and adds a new q-gram
  * right of the block as the window moves along right one position.
  * The best match in the small search space is obtained by running
  * incremental qgram checks to get the best alignment with the lowest error.
  * 
  * @param oldBest The structure holding the best result discovered using
  * the hueristic technique
  * @param XX The qgram vector of xx
  * @param Y The qgram vector of y
  * @return
  */
struct BestMatch hCSC::refine(struct BestMatch oldBest, vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y)
{
    unsigned int i, j, k, h = 0;
    unsigned int startPos, tempStartPos;
    unsigned int startBlock, tempStartBlock, middleBlockIndex;
    unsigned int score, tempScore;
    unsigned int bestScore = UINT_MAX;
    unsigned int bestTPos = oldBest.pos;
    unordered_map<WORD, unsigned int> middleBlock;
    bool jump = false;

    unordered_map<WORD, int> x;
    x.reserve(this->qyNum + this->qxNum);

    unsigned int qgStartPos, qgEndPos;
    string qgStartGram, qgEndGram;
    WORD qgStartWord, qgEndWord;

    if (oldBest.index == 0) {
	startBlock = 0; 
	startPos = 0;
	jump = true;
    } else {
	startPos = oldBest.pos - this->bxSize + 1;
	startBlock = oldBest.index - 1;
    }
    
    if (jump) {
	//find last block on end of x (middle of xx) and save an unaltered copy of it
	middleBlockIndex = (unsigned int)((0.5 * this->m) / (double)this->bxSize);
	middleBlock = XX[middleBlockIndex];
    }

    do {
	//remove qgram left of the block and add q-grams right of the block with every window frame progression @todo Check if this can be refined further
	if (startPos > 0) {
	    tempStartBlock = startBlock;
	    tempStartPos = startPos - 1;

	    for (i = 0; i < this->byNum; i++) {

		qgEndPos = tempStartPos + this->bxSize;
		if (qgEndPos > this->m) {
		    break;
		} else {
		    qgStartPos = tempStartPos;
		    qgStartGram = this->xx.substr(qgStartPos, this->qSize);
		    qgStartWord = this->shuffleOntoWord(qgStartGram);
		    XX[tempStartBlock][qgStartWord]--;
		    if (XX[tempStartBlock][qgStartWord] == 0) {
			XX[tempStartBlock].erase(qgStartWord);
		    }

		    qgEndGram = this->xx.substr(qgEndPos, this->qSize);
		    qgEndWord = this->shuffleOntoWord(qgEndGram);
		    XX[tempStartBlock][qgEndWord]++;
		}

		tempStartPos += this->bxSize;
		tempStartBlock++;
	    }
	}

	//find best score
	score = 0;
	for (i = startBlock, j = 0; i < startBlock + this->byNum; i++, j++) {
	    tempScore = 0;
	    
	    //get qgrams scores for Y block
	    for (auto py = Y[j].begin(); py != Y[j].end(); py++) {
		tempScore += abs(this->findInMap(XX[i], py->first) - this->findInMap(Y[j], py->first));
		x[py->first] = -1;
	    }

	    //get qgrams scores for XX block
	    for (auto px = XX[i].begin(); px != XX[i].end(); px++) {
		if (this->findInMap(x, px->first) != -1) {
		    tempScore += abs(this->findInMap(XX[i], px->first) - this->findInMap(Y[j], px->first));
		}
	    }

	    score += tempScore;

	    x.clear();
	}

	//record best score
	if (score < bestScore) {
	    bestScore = score;
	    bestTPos = startPos;
	}

	//move window along one and update block number
	startPos++;
	h++;

	if (jump && h >= this->bxSize) {
	    startBlock = middleBlockIndex;
	    startPos = startBlock * this->bxSize;
	    XX[middleBlockIndex] = middleBlock;
	    jump = false;
	}
    }
    while (h < (2 * this->bxSize - 2) && (startPos + this->n) <= this->m);

    struct BestMatch best;
    best.index = (int)(bestTPos / this->bxSize); //block index starting from 0
    best.score = bestScore; //score
    best.pos = bestTPos; //index

    return best;
}

////////////////////
// Public Methods //
////////////////////

/**
  * Execute heuristic algorithm
  * 
  * @param best The best match
  */
int hCSC::run(struct BestMatch * best)
{
    //check block and q length
    if (this->qSize >= min(this->bxSize, this->bySize)) {
	cerr << "An error occured. Q-gram length cannot be larger or equal to block length." << endl;
	return EXIT_FAILURE;
    }

    //check there isn't any q-gram overflow
    if ((this->qSize * this->l) > WORD_SIZE) {
	cerr << "An error occured. Q-gram overflow exception. Please choose a smaller q-gram size." << endl;
	return EXIT_FAILURE;
    }

    //create XX and Y matrices
    vector<unordered_map<WORD, unsigned int>> XX;
    vector<unordered_map<WORD, unsigned int>> Y;
    XX.resize(this->bxNum);
    Y.resize(this->byNum);
    unsigned int i;
    for (i = 0; i < this->bxNum; i++) {
	XX[i].reserve(this->qxNum);
    }
    for (i = 0; i < this->byNum; i++) {
	Y[i].reserve(this->qyNum);
    }

    XX = fillQGramBlocks(XX, this->xx, this->m, this->bxSize, this->qxNum);
    Y = fillQGramBlocks(Y, this->y, this->n, this->bySize, this->qyNum);

    //compare XX q-gram signature against Y q-grams to find best match
    struct BestMatch bestSoFar = this->getBestScoringBlock(XX, Y);
    //cout << "Heuristic best score, index and position " << bestSoFar.score << " " << bestSoFar.index << " " << bestSoFar.pos << endl;

    //get refinement
    ( * best ) = this->refine(bestSoFar, XX, Y);
    //cout << "Refined best position and score - " << best.pos << ", " << best.score << " - " << this->xx.substr(best.pos, this->n) << endl;

    return EXIT_SUCCESS;
}
