/******************************************************************************
 *  Copyright (C) 2021        Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef BIT_PARALLEL_ED_H
#define BIT_PARALLEL_ED_H

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// ============================================================================
// DEFINITIONS
// ============================================================================

#define WORD_SIZE (64ull)
#define BLOCK_SIZE (32ull)
#define MAX_ED ((WORD_SIZE - BLOCK_SIZE - 2ull) / 3ull)
#define LEFT (2ull * MAX_ED + 1ull)
#define DIAG_R0 (2ull * MAX_ED)

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

typedef struct {
    uint64_t HP;     // bit vector to indicate which delta_H == +1
    uint64_t HN;     // bit vector to indicate which delta_H == -1
    uint64_t VP;     // bit vector to indicate which delta_V == +1
    uint64_t VN;     // bit vector to indicate which delta_V == -1
    uint64_t D0;     // bit vector to indicate which delta_D == 0
    uint64_t RAC;    // bit vector to indicate the Rightmost Active Column
                     // = rightmost column with a value <= maxED
    uint64_t score;  // minScore at the diagonal
} BitVectors;

typedef struct {
    uint64_t A;  // match vector to indicate occurrences of A in X
    uint64_t C;  // match vector to indicate occurrences of C in X
    uint64_t G;  // match vector to indicate occurrences of G in X
    uint64_t T;  // match vector to indicate occurrences of T in X
} MatchVectors;

// ============================================================================
// CLASS BIT-PARALLEL-ED MATRIX
// ============================================================================

class BitParallelED {
   public:
    /**
     * Constructor
     */
    BitParallelED() {
        // create the alphabet mapping
        char2idx = std::vector<char>(256, 4);
        char2idx['A'] = 0;
        char2idx['C'] = 1;
        char2idx['G'] = 2;
        char2idx['T'] = 3;
    }

    /**
     * Bit-encode the horizontal sequence X. Call this routine BEFORE calling
     * initializeMatrix(). You may call initializeMatrix() multiple times
     * (with different initialization settings) with a fixed sequence X.
     */
    void setSequence(const std::string& X) {
        n = X.size() + 1;    // number of columns
        m = 2 * MAX_ED + n;  // this is an upper bound, the exact maxED
                             // is specified during initializeMatrix()

        // allocate and initialize the match vectors
        mv.resize((m + BLOCK_SIZE - 1) / BLOCK_SIZE);

        // encode the first block
        const uint64_t init = (1ull << LEFT) - 1;
        mv[0] = {init, init, init, init};
        uint64_t bitmask = 1ull << LEFT;
        size_t je = std::min<size_t>(X.size(), WORD_SIZE - LEFT);
        for (size_t j = 0; j < je; j++) {
            assert(char2idx[X[j]] < 4);  // assert ACTG alphabet
            mv[0][char2idx[X[j]]] |= bitmask;
            bitmask <<= 1;
        }

        // encode the remaining blocks
        for (size_t b = 1; b < mv.size(); b++) {
            mv[b][0] = mv[b - 1][0] >> BLOCK_SIZE;
            mv[b][1] = mv[b - 1][1] >> BLOCK_SIZE;
            mv[b][2] = mv[b - 1][2] >> BLOCK_SIZE;
            mv[b][3] = mv[b - 1][3] >> BLOCK_SIZE;

            bitmask = 1ull << (WORD_SIZE - BLOCK_SIZE);
            size_t jb = WORD_SIZE - LEFT + (b - 1) * BLOCK_SIZE;
            size_t je = std::min<size_t>(X.size(), jb + BLOCK_SIZE);
            for (size_t j = jb; j < je; j++) {
                assert(char2idx[X[j]] < 4);  // assert ACTG alphabet
                mv[b][char2idx[X[j]]] |= bitmask;
                bitmask <<= 1;
            }
        }
    }

    /**
     * Initialize the alignment matrix
     * @param maxED Maximum edit distance allowed during alignment
     * @param initED Edit distances of column zero (default = 0, 1, ... maxED)
     */
    void initializeMatrix(uint maxED, const std::vector<uint>& initED = {}) {
        // make sure maxED is within supported range
        assert(maxED <= MAX_ED);

        // sanity check on the initED vector
        assert(initED.empty() || initED.front() <= maxED);
        assert(initED.empty() || initED.back() <= maxED);

        this->maxED = maxED;             // store the maximum ED
        Wv = (initED.empty()) ? maxED :  // vertical width of the band
                 initED.size() - 1 + maxED - initED.back();

        // sanity check on the size of initED
        assert(Wv <= 2 * MAX_ED);

        m = Wv + n;  // number of rows

        // allocate and initialize bit vectors
        bv.resize(m);
        bv[0].score = initED.empty() ? 0 : initED[0];
        Wh = maxED - bv[0].score;  // horizontal width of the band

        // initialize top row as [2*MAX_ED, ..., 2, 1, 0, 1, 2, ...]
        bv[0].HP = (~0ull) << LEFT;
        bv[0].HN = ~bv[0].HP;

        // correct top row if initED has been specified
        for (size_t i = 1; i < std::min<size_t>(initED.size(), LEFT + 1); i++) {
            if (initED[i] < initED[i - 1]) {
                bv[0].HP ^= 1ull << (LEFT - i);  // set HP to 1
                bv[0].HN ^= 1ull << (LEFT - i);  // set HN to 0
            } else if (initED[i] == initED[i - 1]) {
                bv[0].HN ^= 1ull << (LEFT - i);  // set HN to 0
            }
        }

        // RAC equals the right-most active element
        bv[0].RAC = 1ull << (DIAG_R0 + Wh);
    }

    /**
     * Initialize the alignment matrix for in-text validation
     * @param maxED Maximum edit distance allowed during alignment
     */
    void initInTextValidation(uint maxED) {
        // make sure maxED is within supported range
        assert(maxED <= MAX_ED);

        this->maxED = maxED;  // store the maximum ED
        Wv = maxED;           // vertical width of the band
        m = Wv + n;           // number of rows

        // allocate and initialize bit vectors
        bv.resize(m);
        bv[0].score = 0;
        Wh = 2 * maxED;  // horizontal width of the band

        // initialize top row as [2*MAX_ED, ..., 2, 1, 0, ..., 0, 1, 2,...]
        bv[0].HP = (~0ull) << (LEFT + maxED);
        bv[0].HN = ~(bv[0].HP >> maxED);

        // RAC equals the right-most active element
        bv[0].RAC = 1ull << (DIAG_R0 + Wh);
    }

    /**
     * Compute a row of the edit distance matrix in a bit-parallel manner
     * @param i row index in range [1...m[
     * @param Y character of Y-sequence at row i
     * @return false if all elements on row i exceed maxED, true otherwise
     */
    bool computeRow(uint i, char Y) {
        assert(i > 0);
        assert(i < m);
        assert(char2idx[Y] < 4);

        // define BLOCK_SIZE as power of two to make sure this is fast:
        const uint b = i / BLOCK_SIZE;  // block identifier
        const uint l = i % BLOCK_SIZE;  // leftmost relevant bit

        // aliases to the bit vectors of the current row i (will be computed)
        uint64_t& HP = bv[i].HP;
        uint64_t& HN = bv[i].HN;
        uint64_t& VP = bv[i].VP;
        uint64_t& VN = bv[i].VN;
        uint64_t& D0 = bv[i].D0;
        uint64_t& RAC = bv[i].RAC;

        // select the right match vector
        const uint64_t& M = mv[b][char2idx[Y]];

        // copy the input vectors pertaining the previous row i-1
        HP = bv[i - 1].HP;
        HN = bv[i - 1].HN;
        RAC = bv[i - 1].RAC << 1;

        // if we are entering a new block, shift input vectors to the right
        // so that they align with the current block
        if (i % BLOCK_SIZE == 0) {
            HP >>= BLOCK_SIZE;
            HN >>= BLOCK_SIZE;
            RAC >>= BLOCK_SIZE;
        }

        // compute the 5 bitvectors that encode the edit distance minScore (Hyyro)
        D0 = (((M & HP) + HP) ^ HP) | M | HN;
        VP = HN | ~(D0 | HP);
        VN = D0 & HP;
        HP = (VN << 1) | ~(D0 | (VP << 1));
        HN = (D0 & (VP << 1));

        // compute the minScore at the diagonal
        const size_t diagBit = l + DIAG_R0;
        bv[i].score = bv[i - 1].score + (D0 & (1ull << diagBit) ? 0 : 1);

        // update the rightmost active column (Hyyro)
        if ((D0 & RAC) == 0) {
            size_t val = 1u;
            while (val > 0) {
                if (HP & RAC) val--;
                if (HN & RAC) val++;
                if (RAC == (1ull << (diagBit - Wv))) return false;
                RAC >>= 1;
            }
        }

        return true;
    }

    /**
     * Find the minimum edit distance value and its position on a row
     * @param i Row index
     * @param jMin Column index at which minimum value is found (output)
     * @param minScore Mimumum value (output)
     */
    void findMinimumAtRow(uint i, uint& jMin, uint& minScore) const {
        // assume the minimal value is found at the diagonal
        minScore = bv[i].score;
        jMin = i;

        // check for lower values left of the diagonal
        uint thisScore = bv[i].score;
        uint64_t bit = 1ull << ((i % BLOCK_SIZE) + DIAG_R0);
        for (uint j = i; j >= getFirstColumn(i); j--, bit >>= 1) {
            if (bv[i].HP & bit) thisScore--;
            if (bv[i].HN & bit) thisScore++;
            if (thisScore < minScore) {
                minScore = thisScore;
                jMin = j - 1;
            }
        }

        // check values to the right of the diagonal
        thisScore = bv[i].score;
        bit = 1ull << ((i % BLOCK_SIZE) + LEFT);
        for (uint j = i + 1; j <= getLastColumn(i); j++, bit <<= 1) {
            if (bv[i].HP & bit) thisScore++;
            if (bv[i].HN & bit) thisScore--;
            if (thisScore < minScore) {
                minScore = thisScore;
                jMin = j;
            }
        }
    }

    /**
     * Do backtracking and compute CIGAR string
     * @param Q Query seq.  Ref. seq. should be set using setSequenceX(...)
     * @param refBegin Begin offset of the reference sequence (output)
     * @param refEnd End offset of the reference sequence (output)
     * @param ED Edit distance score associated with this alignment (output)
     * @param CIGAR CIGAR string (output)
     */
    void trackBack(const std::string& Q, uint& refBegin, uint& refEnd, uint& ED,
                   std::vector<std::pair<char, uint>>& CIGAR) const {
        // sanity check
        assert(Q.size() < m);

        CIGAR.clear();
        CIGAR.reserve(2 * MAX_ED + 1);

        // find the mimimum ED score at the end of the query sequence
        uint i = (uint)Q.size();
        findMinimumAtRow(i, refEnd, ED);

        // backtracking
        uint j = refEnd;
        while (i > 0) {
            const uint b = i / BLOCK_SIZE;  // block identifier
            const uint64_t& M = mv[b][char2idx[Q[i - 1]]];
            uint64_t bit = 1ull << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if (j && ((M | ~(M | bv[i].D0)) & bit)) {  // diagonal
                i--;
                j--;
                if (CIGAR.empty() || CIGAR.back().first != 'M')
                    CIGAR.push_back(std::make_pair('M', 1));
                else
                    CIGAR.back().second++;
            } else if (j && (bv[i].HP & bit)) {  // horizontal gap
                j--;
                if (CIGAR.empty() || CIGAR.back().first != 'D')
                    CIGAR.push_back(std::make_pair('D', 1));
                else
                    CIGAR.back().second++;
            } else if (bv[i].VP & bit) {  // vertical gap
                i--;
                if (CIGAR.empty() || CIGAR.back().first != 'I')
                    CIGAR.push_back(std::make_pair('I', 1));
                else
                    CIGAR.back().second++;
            }
        }

        std::reverse(CIGAR.begin(), CIGAR.end());
        refBegin = j;
    }

    /**
     * Operator () overloading -- this procedure is O(|i-j|)
     * @param i Row index
     * @param j Column index
     * @return Score at position (i, j)
     */
    uint operator()(uint i, uint j) const {
        assert(i < m);  // make sure i and j are within matrix bounds
        assert(j < n);

        uint score = bv[i].score;
        const uint bit = (i % BLOCK_SIZE) + DIAG_R0;
        if (i > j) {
            for (uint o = 0; o < i - j; o++) {
                score -= (bv[i].HP >> (bit - o)) & 1ull;
                score += (bv[i].HN >> (bit - o)) & 1ull;
            }
        } else {
            for (uint o = 1; o <= j - i; o++) {
                score += (bv[i].HP >> (bit + o)) & 1ull;
                score -= (bv[i].HN >> (bit + o)) & 1ull;
            }
        }

        return score;
    }

    /**
     * Check whether after row i, the alignment involves only vertical gaps.
     * This happens when row i includes the final column n and when all values
     * on row i decrease monotically
     * @return true of false
     */
    bool onlyVerticalGapsLeft(uint i) const {
        assert(i < m);

        if (i + LEFT < n)  // if the column n is not yet reached on row i
            return false;

        const uint b = i / BLOCK_SIZE;
        const uint r = i % BLOCK_SIZE;

        // check if all relevant bits for HN are set to 1
        size_t bb = DIAG_R0 - Wv + r + 1;
        size_t be = DIAG_R0 + n - b * BLOCK_SIZE;
        return (((~bv[i].HN >> bb) << bb) << (64 - be)) == 0ull;
    }

    /**
     * Retrieves the first column index that needs to be filled in for the row
     * @param i the row to fill in
     * @returns the first column to fill in
     */
    uint getFirstColumn(uint i) const { return (i <= Wv) ? 1u : i - Wv; }

    /**
     * Retrieves the last column index that needs to be filled in for the row
     * @param i the row to fill in
     * @returns the last column to fill in
     */
    uint getLastColumn(uint i) const { return std::min(n - 1, i + Wh); }

    /**
     * Get the number of columns in the matrix
     * @return The number of columns in the matrix (== X.size() + 1)
     */
    uint getNumberOfCols() const { return n; }

    /**
     * Get the number of rows in the matrix
     * @return The number of rows in the matrix (== Y.size() + 1)
     */
    uint getNumberOfRows() const { return m; }

    /**
     * Check whether setSequenceX() has been called
     * @return True or false
     */
    bool sequenceSet() const { return !mv.empty(); }

    /**
     * Get the vertical size of the final column
     * @return True or false
     */
    uint getSizeOfFinalColumn() const { return Wh + Wv + 1; }

    /**
     * Print the banded matrix
     * @param maxRow Last row index to print
     */
    void printMatrix(uint maxRow = 500) const {
        for (uint i = 0; i < std::min<uint>(maxRow + 1, m); i++) {
            uint firstCol = (i <= Wv) ? 0u : getFirstColumn(i);
            uint lastCol = getLastColumn(i);
            std::cout << (i < 10 ? "0" : "") << std::to_string(i);
            std::cout << " [" << getFirstColumn(i) << "," << getLastColumn(i) << "]\t";

            for (uint j = 0; j < firstCol; j++) std::cout << "  ";
            for (uint j = firstCol; j <= lastCol; j++) std::cout << operator()(i, j) << " ";

            std::cout << "\tRAC:" << -(int)Wv + (int)i << "/" << std::log2(bv[i].RAC) - DIAG_R0;
            std::cout << (onlyVerticalGapsLeft(i) ? " - true" : " - false");
            uint minScore, minJ;
            findMinimumAtRow(i, minJ, minScore);
            std::cout << "  Min: " << minScore << "@" << minJ;
            std::cout << std::endl;
        }
    }

   private:
    std::vector<char> char2idx;

    uint maxED;  // maximum allowed edit distance
    uint m;      // number of rows
    uint n;      // number of columns
    uint Wv;     // vertical width of the band
    uint Wh;     // horizontal width of the band

    std::vector<BitVectors> bv;               // bit vectors
    std::vector<std::array<uint64_t, 4>> mv;  // match vectors
};

#endif
