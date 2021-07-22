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

#include "bitpared.h"

#include <fstream>
#include <cstdio>

#include "gtest/gtest.h"

using namespace std;

TEST(BitParEDTest, GlobalAlignment)
{
    // AAACCCGGGTATTT-AT
    //     *
    // AAACAC-GGTATTTAAT
    string X = "AAACCCGGGTATTTAT";
    string Y = "AAACACGGTATTTAAT";

    // basic global alignment demo
    BitParallelED myBPED;
    myBPED.setSequence(X);                      // set and bit-encode sequence X
    myBPED.initializeMatrix(4);                 // initialize matrix with maxED
    for (uint i = 1; i <= Y.size(); i++) {      // compute all rows
        bool valid = myBPED.computeRow(i, Y[i-1]);
        EXPECT_EQ(valid, true);
    }
    uint score = myBPED(Y.size(), X.size());    // get score
    EXPECT_EQ(score, 3);

    // CIGAR string computation
    vector<pair<char, uint> > CIGAR;
    uint refBegin, refEnd;
    myBPED.trackBack(Y, refBegin, refEnd, score, CIGAR);

    EXPECT_EQ(refBegin, 0u);
    EXPECT_EQ(refEnd, X.size());
    EXPECT_EQ(score, 3);
    ASSERT_EQ(CIGAR.size(), 5);
    EXPECT_EQ(CIGAR[0].first, 'M'); EXPECT_EQ(CIGAR[0].second, 6);
    EXPECT_EQ(CIGAR[1].first, 'D'); EXPECT_EQ(CIGAR[1].second, 1);
    EXPECT_EQ(CIGAR[2].first, 'M'); EXPECT_EQ(CIGAR[2].second, 7);
    EXPECT_EQ(CIGAR[3].first, 'I'); EXPECT_EQ(CIGAR[3].second, 1);
    EXPECT_EQ(CIGAR[4].first, 'M'); EXPECT_EQ(CIGAR[4].second, 2);
}
