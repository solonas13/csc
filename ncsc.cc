
#include "ncsc.h"

/////////////////////
// Private Methods //
/////////////////////

/**
  * Run Naive: Run y against x and move y one position along for n positions
  * 
  * @param XX The qgram vector of xx
  * @param Y The qgram vector of y
  * @return The best position of y in xx
  */
struct BestMatch nCSC::runNaive(vector<unordered_map<WORD, unsigned int>> XX, vector<unordered_map<WORD, unsigned int>> Y)
{
    unsigned int i, j;
    unsigned int startPos = 0, startBlock = 0;
    unsigned int tempStartPos, tempStartBlock;
    unsigned int score, tempScore;
    unsigned int bestScore = UINT_MAX, bestTPos;

    unordered_map<WORD, int> x;
    x.reserve(this->qyNum + this->qxNum);

    unsigned int qgStartPos, qgEndPos;
    string qgStartGram, qgEndGram;
    WORD qgStartWord, qgEndWord;

    do {
	//remove qgram left of the block and add q-grams right of the block with every window frame progression
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
    }
    while (startPos < this->n);

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
  * Execute naive algorithm
  * 
  * @param best The best match
  */
int nCSC::run(struct BestMatch * best)
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

    //print out all scores of y against x
    ( * best ) = this->runNaive(XX, Y);
    //cout << "The best best position and score of y in xx is - " << best.pos << ", " << best.score << " - " << this->xx.substr(best.pos, this->n) << endl;

    return EXIT_SUCCESS;
}
