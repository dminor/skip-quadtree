/*
Copyright (c) 2010 Daniel Minor 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>

#include <ANN/ANN.h>

int main(int argc, char **argv)
{ 
    if (argc < 2) {
        std::cout << "usage: knn <pts> [queries] [epsilon]" << std::endl;
        exit(1);
    }

    //read points
    int pt_count;

    std::ifstream ptf(argv[1]);

    if (!ptf) {
        std::cout << "error: could not open points file: " << argv[1] << std::endl;
        exit(1); 
    }

    ptf >> pt_count;

    if (pt_count < 0) {
        std::cout << "error: invalid point count " << pt_count << std::endl;
        exit(1);
    }

    ANNpointArray pts = annAllocPts(pt_count, 2); 
    for (int i = 0; i < pt_count; ++i) { 
        double d;
        char c;

        ptf >> d; pts[i][0] = d;
        ptf >> c;
        ptf >> d; pts[i][1] = d; 
    }

    ptf.close();
    
    ANNkd_tree kt(pts, pt_count, 2);
    if (argc < 3) return 1;

    //read queries
    int q_count;

    std::ifstream qf(argv[2]);

    if (!qf) {
        std::cout << "error: could not open query file: " << argv[2] << std::endl;
        exit(1); 
    }

    qf >> q_count;

    if (q_count < 0) {
        std::cout << "error: invalid query count " << q_count << std::endl;
        exit(1);
    }

    ANNpointArray queries = annAllocPts(q_count, 2);

    for (int i = 0; i < q_count; ++i) {
        double d;
        char c;

        qf >> d; queries[i][0] = d;
        qf >> c;
        qf >> d; queries[i][1] = d;
    }

    qf.close();

    //read query epsilon
    double epsilon = 0.0;
    if (argc == 4) epsilon = atof(argv[3]);

    //run queries
    ANNidx *nn_idx = new ANNidx[5];
    ANNdist *dists = new ANNdist[5];

    for (int i = 0; i < q_count; ++i) { 

        kt.annkSearch(queries[i], 5, nn_idx, dists, epsilon);

        std::cout << "query " << i << ": (" << queries[i][0] << ", " << queries[i][1] << ")" << std::endl;

        for (int j = 0; j < 5; ++j) { 
            std::cout << "(" << pts[nn_idx[j]][0] << ", " << pts[nn_idx[j]][1] << ") " << dists[j] << std::endl; 
        } 
    }

    std::cout << "done." << std::endl;

    delete[] nn_idx;
    delete[] dists;
    annClose();

    return 0;
}
