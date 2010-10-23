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

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "compressed_quadtree.h"

struct Point {
    double coords[2];

    Point() {};

    Point(const Point &other)
    {
        coords[0] = other.coords[0];
        coords[1] = other.coords[1]; 
    }

    double operator[](size_t idx) const {return coords[idx];}
    double &operator[](size_t idx) {return coords[idx];}
};

double query_x;
double query_y;

bool pt_lt(const Point pt1, const Point pt2)
{ 
    double d1 = (pt1[0]-query_x)*(pt1[0]-query_x)+(pt1[1]-query_y)*(pt1[1]-query_y);
    double d2 = (pt2[0]-query_x)*(pt2[0]-query_x)+(pt2[1]-query_y)*(pt2[1]-query_y);

    return d1 < d2;
}

std::vector<Point *> sort_knn(Point *pts, size_t pt_count, size_t k, Point pt)
{
    std::vector<Point *> qr;

    query_x = pt[0]; query_y = pt[1];
    std::sort(pts, pts + pt_count, pt_lt);

    for (int i = 0; i < k; ++i) {
        qr.push_back(&pts[i]);
    } 

    return qr;
}

int main(int argc, char **argv)
{ 
    if (argc != 3) {
        std::cout << "usage: knn <pts> <queries>" << std::endl;
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

    Point *pts = new Point[pt_count]; 
    for (int i = 0; i < pt_count; ++i) { 
        char c;
        ptf >> (pts[i][0]);
        ptf >> c;
        ptf >> (pts[i][1]);
    }

    ptf.close();

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

    Point *queries = new Point[q_count]; 

    for (int i = 0; i < q_count; ++i) {
        char c;
        qf >> queries[i][0];
        qf >> c;
        qf >> queries[i][1];
    }

    qf.close();

    CompressedQuadtree<Point> cqt(2, pts, pt_count);

    //run queries
    for (int i = 0; i < q_count; ++i) { 

        std::vector<std::pair<Point *, double> > sqr = cqt.knn(5, queries[i], 0.0);  
        std::vector<Point *> lqr = sort_knn(pts, pt_count, 5, queries[i]);  

        for (int j = 0; j < 5; ++j) {

            if (sqr[j].first == 0) {
                std::cout << "error: did not find 5 neighbouring points in tree" << std::endl;
                return 1; 
            }

            double x = (*lqr[j])[0] - (*sqr[j].first)[0]; 
            double y = (*lqr[j])[1] - (*sqr[j].first)[1];

            if ((x*x + y*y) > 0.0001) {
                std::cout << "error: compressed quadtree nearest neighbour does not match sort nearest neighbour" << std::endl;

                std::cout << "query " << i << ": " << queries[i][0] << " " << queries[i][1] << std::endl;
                double d2 = (queries[i][0]-(*lqr[j])[0])*(queries[i][0]-(*lqr[j])[0]) + (queries[i][1]-(*lqr[j])[1])*(queries[i][1]-(*lqr[j])[1]);
                std::cout << "result: " << j << std::endl;
                std::cout << "sqr: " <<  (*sqr[j].first)[0] << " " << (*sqr[j].first)[1] << " " << sqr[j].second << std::endl;
                std::cout << "lqr: " << (*lqr[j])[0] << " " << (*lqr[j])[1] << " " << sqrt(d2) << std::endl;

                return 1;
            }

        } 
    }

    std::cout << "done." << std::endl;

    delete[] pts;
    delete[] queries; 

    return 0;
}
