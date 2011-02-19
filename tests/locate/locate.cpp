/*
Copyright (c) 2011 Daniel Minor 

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
#include "skip_quadtree.h"

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

int main(int argc, char **argv)
{ 
    if (argc < 2) {
        fprintf(stderr, "usage: locate <pts> \n");
        exit(1);
    }

    //read points
    int pt_count;

    FILE *f = fopen(argv[1], "r");

    if (!f) {
        printf("error: could not open points file: %s\n", argv[1]);
        return -1;
    }

    int result = fscanf(f, "%d", &pt_count);

    if (result != 1|| pt_count < 0) {
        fprintf(stderr, "error: invalid point count %d\n", pt_count);
        return -1;
    }

    Point *pts = new Point[pt_count];

    double x, y;
    for (int i = 0; i < pt_count; ++i) {
        int result = fscanf(f, "%lf, %lf", &x, &y);
        if (result != 2) {
            fprintf(stderr, "error: could not read point %d\n", i);
            delete[] pts;
            return -1;
        }

        pts[i][0] = x;
        pts[i][1] = y; 
    }

    fclose(f);

    CompressedQuadtree<Point> cqt(2, pts, pt_count);

    for (size_t i = 0; i < pt_count; ++i) {
        CompressedQuadtree<Point>::Node *node = cqt.locate(pts[i]);
        if (node) {
            for (size_t d = 0; d < 2; ++d) { 
                if (node->mid[d] - node->side[d] - pts[i][d] > 0.1 || pts[i][d] - node->mid[d] - node->side[d] > 0.1) { 
                    fprintf(stderr, "error: compressed quadtree found node not containing point %d: (%.1lf, %.1lf)\n", i, pts[i][0], pts[i][1]);
                    return -1;
                } 
            } 
        } else {
            fprintf(stderr, "error: compressed quadtree failed to locate point %d: (%.1lf, %.1lf)\n", i, pts[i][0], pts[i][1]);
            return -1;
        }
    }

    SkipQuadtree<Point> sqt(2, pts, pt_count);
    for (size_t i = 0; i < pt_count; ++i) {
        SkipQuadtree<Point>::Node *node = sqt.locate(pts[i]);
        if (node) {
            for (size_t d = 0; d < 2; ++d) { 
                if (node->mid[d] - node->side[d] - pts[i][d] > 0.1 || pts[i][d] - node->mid[d] - node->side[d] > 0.1) { 
                    fprintf(stderr, "error: skip quadtree found node not containing point %d: (%.1lf, %.1lf)\n", i, pts[i][0], pts[i][1]);
                    return -1;
                } 
            } 
        } else {
            fprintf(stderr, "error: skip quadtree failed to locate point %d: (%.1lf, %.1lf)\n", i, pts[i][0], pts[i][1]);
            return -1;
        } 
    }

    delete[] pts;
}
