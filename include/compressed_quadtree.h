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

#ifndef COMPRESSED_QUADTREE_H_
#define COMPRESSED_QUADTREE_H_

/*
    Compressed Quadtree implementation supporting approximate nearest neighbour
    queries, based upon the description in:  
    
    Eppstein, D., Goodrich, M. T., Sun, J. Z. (2008) The Skip Quadtree:
    A Simple Dynamic Data Structure for Multidimensional Data,
    Int. Journal on Computational Geometry and Applications, 18(1/2), pp. 131 - 160
*/

#include <algorithm>
#include <limits>
#include <vector>
#include <iostream>
#include <cstring>

template<class Point> class CompressedQuadtree {

    public:

        //Nodes of the quadtree
        struct Node { 
            Node **nodes;       //children
            Point mid;          //midpoint
            Point side;         //half side length
            Point *pt;          //point, if data stored
        };

        //Keep track of point and priority for nearest neighbour search
        struct NodeDistance {

            Node *node;
            double distance;

            NodeDistance(Node *node, double distance) : node(node), distance(distance) {};

            //min-heap
            bool operator<(const NodeDistance &other) const
            {
                return distance > other.distance;
            }
        }; 

        CompressedQuadtree(size_t dim, Point *pts, size_t n) : dim(dim), nnodes(1 << dim)
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

        virtual ~CompressedQuadtree()
        {
           if (root) delete_worker(root);  
        }

        std::vector<std::pair<Point *, double> > knn(size_t k, const Point &pt, double eps) 
        {
            //setup query result vector
            std::vector<std::pair<Point *, double> > qr; 
            qr.reserve(k);
            for (size_t i = 0; i < k; ++i) {
                qr.push_back(std::make_pair<Point *, double>(0, std::numeric_limits<double>::max()));
            }
 
            //initialize priority queue for search
            std::vector<NodeDistance > pq; 
            pq.push_back(NodeDistance(root, 0.0));

            while (!pq.empty()) {

                std::pop_heap(pq.begin(), pq.end());
                Node *node = pq.back().node;
                double node_dist = pq.back().distance; 
                pq.pop_back();

                if (node->nodes == 0) {

                    //calculate distance from query point to this point
                    double dist = 0.0; 
                    for (size_t d = 0; d < dim; ++d) {
                        dist += ((*node->pt)[d]-pt[d]) * ((*node->pt)[d]-pt[d]); 
                    }

                    //insert point in proper order, for small k this will be fast enough
                    for (size_t i = 0; i < k; ++i) {
                        if (qr[i].second > dist) {
                            for (size_t j = k - 1; j > i; --j) {
                                qr[j] = qr[j - 1];
                            }
                            qr[i].first = node->pt;
                            qr[i].second = dist;
                            break; 
                        }
                    } 
                } else {

                    //find k-th distance
                    double kth_dist = qr[k - 1].second; 

                    //stop searching, all further nodes farther away than k-th value
                    if (kth_dist <= (1.0 + eps)*node_dist) {
                        break;
                    }

                    for (size_t n = 0; n < nnodes; ++n) { 
                        //calculate distance to each of the non-zero children
                        //if less than k-th distance, then visit 
                        if (node->nodes[n]) {

                            double min_dist = min_pt_dist_to_node(pt, node->nodes[n]);

                            //if closer than k-th distance, search
                            if (min_dist < kth_dist) { 
                                pq.push_back(NodeDistance(node->nodes[n], min_dist));
                                std::push_heap(pq.begin(), pq.end()); 
                            }
                        } 
                    }
                }
            } 

            return qr;
        }

        Node *locate(const Point &pt) 
        {

            Node *node = 0;

            //see if within bounds covered by the quadtree 
            bool in_quadtree = true;

            for (size_t d = 0; d < dim; ++d) { 
                if (root->mid[d] - root->side[d] >= pt[d] || pt[d] >= root->mid[d] + root->side[d]) { 
                    in_quadtree = false; 
                    break; 
                } 
            }

            //search for node containing the query point 
            if (in_quadtree) { 
                node = root; 

                while (node) { 
                    if (node->nodes) { 
                        size_t n = 0; 
                        for (size_t d = 0; d < dim; ++d) { 
                            if (pt[d] > node->mid[d]) n += 1 << d; 
                        } 

                        if (node->nodes[n]) node = node->nodes[n]; 
                        else break;

                    } else { 
                        break; 
                    } 
                } 
            } 

            return node; 
        }

        Node *root;

    private:

        size_t dim; 
        size_t nnodes;

        Node *worker(const Point &mid, const Point &side, std::vector<Point *> &pts)
        {
            Node *node = new Node; 
            for (size_t d = 0; d < dim; ++d) {
                node->mid[d] = mid[d];
                node->side[d] = side[d]; 
            }

            if (pts.size() == 1) {
                node->nodes = 0;
                node->pt = pts[0];
            } else { 
                node->nodes = new Node *[nnodes];
                node->pt = 0;

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
                            Node *temp = node;
                            node = node->nodes[n];
                            delete temp;
                            break;
                        }
                    }
                }
            } 

            return node;
        } 

        double min_pt_dist_to_node(const Point &pt, Node *node)
        {
            bool inside = true; 
            double min_dist = std::numeric_limits<double>::max();
            for (size_t d = 0; d < dim; ++d) { 
        
                double dist; 
                if (pt[d] < node->mid[d] - node->side[d]) { 
                    dist = node->mid[d] - node->side[d] - pt[d];
                    inside = false;
                } else if (pt[d] > node->mid[d] + node->side[d]) {
                    dist = pt[d] - node->mid[d] + node->side[d]; 
                    inside = false;
                }

                if (dist < min_dist) min_dist = dist; 
            } 

            if (inside) return 0.0;
            else return min_dist*min_dist;
        }

        void delete_worker(Node *node)
        {
            if (node->nodes) {
                for (size_t n = 0; n < nnodes; ++n) { 
                    if (node->nodes[n]) delete_worker(node->nodes[n]); 
                }

                delete node->nodes;
            }

            delete node; 
        }    

};

#endif

