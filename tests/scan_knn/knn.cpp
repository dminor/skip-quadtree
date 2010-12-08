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
#include <limits>
#include <vector>

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

bool qr_lt(const std::pair<Point *, double> &pt1, const std::pair<Point *, double> &pt2)
{ 
    return pt1.second < pt2.second;
} 

double pt_dist(const Point pt, const Point pt1)
{ 
    return (pt1[0]-pt[0])*(pt1[0]-pt[0])+(pt1[1]-pt[1])*(pt1[1]-pt[1]);
}

//
// Return k nearest neighbours by linearly scanning input
//
std::vector<std::pair<Point *, double> > scan_knn(Point *pts, size_t pt_count, size_t k, Point pt)
{
    std::vector<std::pair<Point *, double> > qr;

    double kth_dist = std::numeric_limits<double>::min();
    for (size_t i = 0; i < pt_count; ++i) {
        double d = pt_dist(pt, pts[i]);

        if (qr.size() < k) {
            qr.push_back(std::pair<Point *, double>(&pts[i], d));
            if (d > kth_dist) kth_dist = d;
        } else {
            if (d < kth_dist) {
                double new_kth_dist = std::numeric_limits<double>::min();
                for (size_t j = 0; j < k; ++j) {
                    if (qr[j].second == kth_dist) {
                        qr[j].first = &pts[i];
                        qr[j].second = d;
                    }

                    if (qr[j].second > new_kth_dist) new_kth_dist = qr[j].second;
                }

                kth_dist = new_kth_dist;
            }
        }
    }

    //sort query results
    std::sort(qr.begin(), qr.end(), qr_lt);

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

    //run queries
    for (int i = 0; i < q_count; ++i) { 

        std::vector<std::pair<Point *, double> > qr = scan_knn(pts, pt_count, 5, queries[i]);  
        std::cout << "query " << i << ": (" << queries[i][0] << ", " << queries[i][1] << ")" << std::endl;

        for (int j = 0; j < 5; ++j) { 
            std::cout << "(" << (*qr[j].first)[0] << ", " << (*qr[j].first)[1] << ") " << qr[j].second << std::endl; 
        } 
    }

    std::cout << "done." << std::endl;

    delete[] pts;
    delete[] queries; 

    return 0;
}
