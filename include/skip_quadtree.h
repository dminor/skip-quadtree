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

#ifndef SKIP_QUADTREE_H_
#define SKIP_QUADTREE_H_

#include <iostream>
#include <limits>
#include <vector>

template<class Point> class SkipQuadtree {

    public:

        //Nodes of the quadtree
        template<class T> struct Node { 
            Node **nodes;   //children
            T mid;          //midpoint
            T side;         //half side length
            T pt;           //point, if data stored
        };

        //Keep track of point and priority for nearest neighbour search
        template<class T> struct PointPriority {

            T *pt;
            double priority;

            PointPriority(T *pt, double priority) : pt(pt), priority(priority) {};

            //min-heap
            bool operator<(const PointPriority &other) const
            {
                return priority > other.priority;
            }
        }; 

        SkipQuadtree(size_t dim, Point *pts, size_t n) : dim(dim), nnodes(1 << dim)
        { 
            //calculate bounds
            Point bounds[2];
            for (size_t d = 0; d < dim; ++d) {
                bounds[0][d] = std::numeric_limits<double>::max(); 
                bounds[1][d] = std::numeric_limits<double>::min(); 
            }

            for (size_t d = 0; d < dim; ++d) {
                for (size_t i = 0; i < n; ++i) {
                    if (pts[i][d] < bounds[0][d]) bounds[0][d] = pts[i][d];
                    if (pts[i][d] > bounds[1][d]) bounds[1][d] = pts[i][d];
                }
            }
 
            //calculate mid point and half side length
            Point mid, side; 
            for (size_t d = 0; d < dim; ++d) {
                mid[d] = (bounds[0][d]+bounds[1][d]) / 2;
                side[d] = (bounds[1][d]-bounds[0][d]) / 2;
            } 

            //set up points vector 
            std::vector<Point *> pts_vector;
            for (size_t i = 0; i < n; ++i) {
                pts_vector.push_back(&pts[i]);
            }
 
            root = worker(mid, side, pts_vector);
        }

        virtual ~SkipQuadtree()
        { 
        }

        std::vector<PointPriority<Point> > knn(size_t k, const Point &pt) 
        {
            std::vector<PointPriority<Point> > pq; 
            knn_search(pq, root, k, pt);
            std::sort_heap(pq.begin(), pq.end());
            std::reverse(pq.begin(), pq.end());
            return pq;
        }

        Node<Point> *root;

    private:

        size_t dim; 
        size_t nnodes;

        Node<Point> *worker(const Point &mid, const Point &side, std::vector<Point *> &pts)
        {
            Node<Point> *node = new Node<Point>; 
            for (size_t d = 0; d < dim; ++d) {
                node->mid[d] = mid[d];
                node->side[d] = side[d]; 
            }

            if (pts.size() == 1) {
                node->nodes = 0;
                for (size_t d = 0; d < dim; ++d) {
                    node->pt[d] = (*pts[0])[d];
                }
            } else { 

                node->nodes = new Node<Point> *[nnodes];

                //divide points between the nodes 
                std::vector<Point *> *node_pts = new std::vector<Point *>[nnodes];
                for (typename std::vector<Point *>::iterator itor = pts.begin(); itor != pts.end(); ++itor) {

                    //determine node index based upon which which side of midpoint for each dimension
                    size_t n = 0;
                    for (size_t d = 0; d < dim; ++d) {
                        if ((*(*itor))[d] > mid[d]) n += 1 << d; 
                    } 

                    node_pts[n].push_back(*itor);
                } 

                //create new nodes recursively
                size_t ninteresting = 0;
                for (size_t n = 0; n < nnodes; ++n) {

                    if (node_pts[n].size()) {

                        Point nbounds[2];
                        Point new_mid, new_side;
                        for (size_t d = 0; d < dim; ++d) { 
                            new_side[d] = side[d] / 2;
                            if (n & (1 << d)) {
                                new_mid[d] = mid[d] + new_side[d];
                            } else { 
                                new_mid[d] = mid[d] - new_side[d];
                            }
                        }

                        ++ninteresting;
                        node->nodes[n] = worker(new_mid, new_side, node_pts[n]);
                    } else {
                        node->nodes[n] = 0; 
                    }
                }

                delete[] node_pts;

                //compress if less than 2 interesting nodes
                if (ninteresting < 2) {
                    for (size_t n = 0; n < nnodes; ++n) {
                        if (node->nodes[n]) {
                            Node<Point> *temp = node;
                            node = node->nodes[n];
                            delete temp;
                            break;
                        }
                    }
                }
            } 

            return node;
        }

        void knn_search(std::vector<PointPriority<Point> > &pq, Node<Point> *node, size_t k, const Point &pt)
        { 
            if (node->nodes == 0) {
                double d = 0.0; 
                for (int i = 0; i < dim; ++i) {
                    d += (node->pt[i]-pt[i]) * (node->pt[i]-pt[i]); 
                } 
                pq.push_back(PointPriority<Point>(&node->pt, sqrt(d)));
                std::push_heap(pq.begin(), pq.end()); 
            } else {

                //determine node index based upon which which side of midpoint for each dimension
                size_t n_index = 0;
                for (size_t d = 0; d < dim; ++d) {
                    if (pt[d] > node->mid[d]) n_index += 1 << d; 
                } 

                //search for node containing point
                if (node->nodes[n_index]) knn_search(pq, node->nodes[n_index], k, pt);

                //find k-th priority (FIXME: more efficient?) 
                size_t kth_index = 0;
                double kth_dist = 0.0;
                for (int i = 0; i < pq.size() && i < k; ++i) {
                    if (pq[i].priority > kth_dist){
                        kth_dist = pq[i].priority; 
                        kth_index = i;
                    }
                }

                for (size_t n = 0; n < nnodes; ++n) { 
                    //calculate distance to each of the other non-zero nodes
                    //if less than priority, then visit 
                    //update priority, and check next node
                    if (node->nodes[n] && n != n_index) {

                        //find distance to each side of square 
                        for (size_t d = 0; d < dim; ++d) { 
                      
                            //figure out which side of the midpoint are we on 
                            //and calculate distance
                            double dist;
                            if ((*pq[kth_index].pt)[d] < node->nodes[n]->mid[d]) {
                                dist = node->nodes[n]->mid[d] - node->nodes[n]->side[d] - (*pq[kth_index].pt)[d];
                            } else {
                                dist = (*pq[kth_index].pt)[d] - node->nodes[n]->mid[d] + node->nodes[n]->side[d];
                            } 

                            //if closer than k-th distance, search
                            if (dist < kth_dist) {
                                knn_search(pq, node->nodes[n], k, pt);

                                //then need to update priority
                                for (int i = 0; i < pq.size() && i < k; ++i) {
                                    if (pq[i].priority > kth_dist){
                                        kth_dist = pq[i].priority; 
                                        kth_index = i;
                                    }
                                } 
                                break;
                            }
                        }

                    } 

                }
            }
        } 

};

#endif

