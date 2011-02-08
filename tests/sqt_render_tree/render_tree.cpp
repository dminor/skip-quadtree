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

#include "skip_quadtree.h"

#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

typedef double Point[2];

void render_tree(FILE *f, struct SkipQuadtree<Point>::Node *tree, size_t depth)
{

    double x1 = tree->mid[0] - tree->side[0];
    double x2 = tree->mid[0] + tree->side[0];
    double y1 = tree->mid[1] - tree->side[1];
    double y2 = tree->mid[1] + tree->side[1];

    //square
    if (depth < 1) fprintf(f, "4 setlinewidth\n");
    else if (depth < 2) fprintf(f, "3 setlinewidth\n"); 
    else if (depth < 4) fprintf(f, "2 setlinewidth\n");
    else fprintf(f, "1 setlinewidth\n");

    fprintf(f, "%.0f %.0f %.0f %.0f node-bounds\n", x1, x2, y1, y2); 

    if (!tree->nodes) {
        fprintf(f, "%.0f %.0f draw-point\n", (*tree->pt)[0], (*tree->pt)[1]);
    } else { 
        for (int i = 0; i < 4; ++i) {
            if (tree->nodes[i]) render_tree(f, tree->nodes[i], depth+1);
        } 
    }
}

int main(int argc, char **argv)
{ 
    if (argc != 3) {
        printf("usage: render_tree <input> <output>\n");
        exit(1);
    }

    int pt_count;

    FILE *f = fopen(argv[1], "r");

    if (!f) {
        printf("error: could not open points file: %s\n", argv[1]);
        exit(1); 
    }

    fscanf(f, "%d", &pt_count);

    if (pt_count < 0) {
        printf("error: invalid point count %d\n", pt_count);
        exit(1);
    }

    Point *pts = new Point[pt_count];

    double x, y;
    for (int i = 0; i < pt_count; ++i) {
        fscanf(f, "%lf, %lf", &x, &y);
        pts[i][0] = x;
        pts[i][1] = y; 
    }


    fclose(f);

    SkipQuadtree<Point> sqt(2, pts, pt_count);
    printf("completed building tree...\n ");

    f = fopen(argv[2], "w");
    if (!f) {
        printf("error: could not open output file: %s\n", argv[2]);
        exit(1); 
    }

    fprintf(f, "%%%%!PS-Adobe-2.0\n");
    fprintf(f, "%%%%Pages: %d\n", sqt.levels.size()); 
    fprintf(f, "/page-begin {\n");
    fprintf(f, "gsave\n");
    fprintf(f, "} def\n");

    fprintf(f, "/page-end {\n");
    fprintf(f, "   grestore\n");
    fprintf(f, "   showpage\n");
    fprintf(f, "} def\n");

    //define point function for later
    fprintf(f, "/draw-point {\n");
    fprintf(f, "    /y exch def\n");
    fprintf(f, "    /x exch def\n");
    fprintf(f, "    gsave\n");
    fprintf(f, "    newpath\n");
    fprintf(f, "    0.5 0.5 0.7 setrgbcolor\n");
    fprintf(f, "    x y 2 0 360 arc\n");
    fprintf(f, "    closepath\n");
    fprintf(f, "    fill\n");
    fprintf(f, "    newpath\n");
    fprintf(f, "    0.4 setgray\n");
    fprintf(f, "    x y 2 0 360 arc\n");
    fprintf(f, "    closepath\n");
    fprintf(f, "    stroke\n");
    fprintf(f, "    grestore\n");
    fprintf(f, "} def\n");

    //node bounding box
    fprintf(f, "/node-bounds {\n");
    fprintf(f, "    /y2 exch def\n");
    fprintf(f, "    /y1 exch def\n");
    fprintf(f, "    /x2 exch def\n");
    fprintf(f, "    /x1 exch def\n");
    fprintf(f, "    gsave\n");
    fprintf(f, "    0.7 setgray\n");
    fprintf(f, "    newpath\n");
    fprintf(f, "    x2 y2 moveto\n");
    fprintf(f, "    x1 y2 lineto\n");
    fprintf(f, "    x1 y1 lineto\n");
    fprintf(f, "    x2 y1 lineto\n");
    fprintf(f, "    closepath\n");
    fprintf(f, "    stroke \n");
    fprintf(f, "    grestore\n");
    fprintf(f, "} def\n");

    for (size_t i = 0; i < sqt.levels.size(); ++i) {
        fprintf(f, "%%%%Page: %d\n", i + 1);
        fprintf(f, "page-begin\n"); 
        render_tree(f, sqt.levels[i], 0); 
        fprintf(f, "page-end\n");
    }

    fclose(f);

}
