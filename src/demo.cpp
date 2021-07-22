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

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <bitset>

#include "bitpared.h"

using namespace std;

int main()
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
    for (uint i = 1; i <= Y.size(); i++)        // compute all rows
        if (!myBPED.computeRow(i, Y[i-1]))
            break;
    uint score = myBPED(Y.size(), X.size());    // get score

    // CIGAR string computation
    vector<pair<char, uint> > CIGAR;
    uint refBegin, refEnd;
    myBPED.trackBack(Y, refBegin, refEnd, score, CIGAR);

    cout << "Ref. seq.:\t" << X << "\nQuery seq.:\t" << Y << "\n";
    cout << "Edit distance score: " << score << endl;
    cout << "CIGAR string: ";

    for (uint i = 0; i < CIGAR.size(); i++)
        cout << CIGAR[i].first << CIGAR[i].second;
    cout << endl;

    return EXIT_SUCCESS;
}
