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

/*
    Skip Quadtree implementation supporting approximate nearest neighbour
    queries, based upon the description in:  

    Eppstein, D., Goodrich, M. T., Sun, J. Z. (2008)
    The Skip Quadtree: A Simple Dynamic Data Structure for Multidimensional Data,
    Int. Journal on Computational Geometry and Applications, 18(1/2), pp. 131 - 160 
*/

#include <algorithm>
#include <cstring>
#include <cmath>
#include <limits>
#include <list>
#include <vector>

#include <iostream>

const double EQUIDISTANT_EPSILON = 0.001;

template<class Point> class SkipQuadtree {

    public:

        //Nodes of the quadtree
        struct Node { 
            Node *above;        //pointer to same node in n+1 level of skip structure 
            Node *below;        //pointer to same node in n-1 level of skip structure 
            Node *parent;       //pointer to parent
            Node **ancestors;   //pointers to closest ancestors in each direction
            Node **nodes;       //children
            Point mid;              //midpoint
            double radius;             //half side length
            Point *pt;              //point, if data stored

            size_t seq; 

            Node()
            {
                above = below = parent = 0;
                ancestors = nodes = 0;
                pt = 0;
                seq = 0;
            }
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

        SkipQuadtree(size_t dim, Point *pts, size_t n) : dim(dim), nnodes(1 << dim), locate_eps(0.001)
        { 
            //calculate bounds
            Point bounds[2];
            for (size_t d = 0; d < dim; ++d) {
                bounds[0][d] = std::numeric_limits<double>::max(); 
                bounds[1][d] = -std::numeric_limits<double>::max(); 
            }

            for (size_t d = 0; d < dim; ++d) {
                for (size_t i = 0; i < n; ++i) {
                    if (pts[i][d] < bounds[0][d]) bounds[0][d] = pts[i][d];
                    if (pts[i][d] > bounds[1][d]) bounds[1][d] = pts[i][d];
                }
            }
 
            //calculate mid point and half side length
            Point mid; 
            double radius = 0;
            for (size_t d = 0; d < dim; ++d) {
                mid[d] = (bounds[0][d]+bounds[1][d]) / 2;
                double side = (bounds[1][d]-bounds[0][d]) / 2;
                if (side > radius) radius = side;
            } 

            //set up points vector 
            std::vector<Point *> pts_vector;
            for (size_t i = 0; i < n; ++i) {
                pts_vector.push_back(&pts[i]);
            }
 
            root = worker(mid, radius, pts_vector); 

            //set ancestor pointers
            Node **ancestors = new Node *[2*dim]; 
            memset(ancestors, 0, 2*dim*sizeof(Node *));
            assign_ancestors(root, ancestors);
            delete[] ancestors;
 
            //build skip hierarchy
            levels.push_back(root);
            build_skiphierarchy(root); 

            knn_seq = 0;
        }

        virtual ~SkipQuadtree()
        { 
            for(typename std::vector<Node *>::iterator itor = levels.begin(); itor != levels.end(); ++itor) {
                if (*itor) delete_worker(*itor);
            }
        }

        Node *locate(const Point &pt) 
        { 
            Node *node = 0;

            //search for node containing the query point 
            if (in_node(root, pt)) { 
                node = root; 

                while (node) { 
                    if (node->nodes) { 
                        size_t n = 0; 

                        for (size_t d = 0; d < dim; ++d) { 
                            if (pt[d] > node->mid[d]) n += 1 << d; 
                        } 

                        if (node->nodes[n] && in_node(node->nodes[n], pt)) {
                            node = node->nodes[n]; 
                        } else { 
                            if (node->below) node = node->below; 
                            else break; 
                        } 
                    } else { 
                        if (node->below) node = node->below; 
                        else break; 
                    } 
                } 
            } 

            return node; 
        }

        std::list<std::pair<Point *, double> > knn(size_t k, const Point &pt, double eps) 
        {

            //number of nodes to backtrack during knn search
            size_t backtrack_nodes;
            if (eps > 0.0) {
                backtrack_nodes = ceil(log((1.0 + eps)*sqrt((double)dim)/eps));
            } else {
                backtrack_nodes = -1; 
            }

            ++knn_seq;

            //setup query result vector
            std::list<std::pair<Point *, double> > qr; 
 
            //initialize priority queue for search
            std::vector<NodeDistance > pq; 
            pq.push_back(NodeDistance(root, min_pt_dist_to_node(pt, root)));

            ++root->seq;

            while (!pq.empty()) {

                std::pop_heap(pq.begin(), pq.end());
                Node *p = pq.back().node;
                double p_dist = pq.back().distance; 
                pq.pop_back();

                if (p->nodes == 0) {

                    //calculate distance from query point to this point
                    double dist = 0.0; 
                    for (size_t d = 0; d < dim; ++d) {
                        dist += ((*p->pt)[d]-pt[d]) * ((*p->pt)[d]-pt[d]); 
                    }

                    //insert point in result
                    typename std::list<std::pair<Point *, double> >::iterator itor = qr.begin();
                    while (itor != qr.end() && itor->second < dist) {
                        ++itor;
                    }

                    qr.insert(itor, std::make_pair<Point *, double>(p->pt, dist)); 

                    if (qr.size() > k) qr.pop_back();

                } else {

                    //find k-th distance
                    double kth_dist = qr.size() < k? std::numeric_limits<double>::max() : qr.back().second;

                    //stop searching, all further nodes farther away than k-th value
                    if (kth_dist <= (1.0+eps)*p_dist) break;

                    //If distance to furthest corner is <= (1.0 + eps) * distance to nearest
                    //corner, we can safely insert all children and see if we've found k values
                    //not sure if this helps too much for k-nn
                    //if (max_pt_dist_to_node(pt, p) <= (1.0 + eps)*p_dist) break;

                    Node *q = skip_search_equidistant(pt, p, p_dist);

                    for (size_t n = 0; n < nnodes; ++n) {
                        if (q->ancestors[n] && q->ancestors[n]->seq != knn_seq) { 
                            q->ancestors[n]->seq = knn_seq;

                            double min_dist = min_pt_dist_to_node(pt, q->ancestors[n]);

                            if (min_dist*(1.0+eps) < kth_dist) { 
                                pq.push_back(NodeDistance(q->ancestors[n], min_dist));
                                std::push_heap(pq.begin(), pq.end()); 
                            } 
                        } 
                    }

                    for (size_t i = 0; i < backtrack_nodes; ++i) {

                        if (!q || q == p->parent) break;

                        if (q->nodes) { 
                            for (size_t n = 0; n < nnodes; ++n) { 
                                if (q->nodes[n] && q->nodes[n]->seq != knn_seq) {

                                    q->nodes[n]->seq = knn_seq;

                                    double min_dist = min_pt_dist_to_node(pt, q->nodes[n]);

                                    pq.push_back(NodeDistance(q->nodes[n], min_dist));
                                    std::push_heap(pq.begin(), pq.end()); 
                                }
                            }
                        }

                        q = q->parent;

                    } 
                }
            } 

            return qr;
        }

        Node *root; 
        std::vector<Node *> levels;

        size_t knn_seq;
        double locate_eps;

    private:

        size_t dim; 
        size_t nnodes;

        //
        // Worker function to recursively build compressed quadtree
        // 
        Node *worker(const Point &mid, double radius, std::vector<Point *> &pts)
        {
            Node *node = new Node; 

            for (size_t d = 0; d < dim; ++d) {
                node->mid[d] = mid[d];
            }

            node->radius = radius; 

            if (pts.size() == 1) {
                node->nodes = 0;
                node->pt = pts[0];
            } else { 

                node->nodes = new Node *[nnodes];

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
                Node *last_interesting = 0;
                for (size_t n = 0; n < nnodes; ++n) {

                    if (node_pts[n].size()) {

                        Point nbounds[2];
                        Point new_mid;
                        double new_radius = radius / 2;
                        for (size_t d = 0; d < dim; ++d) { 
                            if (n & (1 << d)) {
                                new_mid[d] = mid[d] + new_radius;
                            } else { 
                                new_mid[d] = mid[d] - new_radius;
                            }
                        }

                        ++ninteresting;
                        node->nodes[n] = worker(new_mid, new_radius, node_pts[n]);
                        node->nodes[n]->parent = node;

                        last_interesting = node->nodes[n];
                    } else {
                        node->nodes[n] = 0; 
                    }
                }

                delete[] node_pts;

                //compress if less than 2 interesting nodes
                if (ninteresting < 2) {
                    delete node;
                    node = last_interesting; 
                }
            } 

            return node;
        }

        //
        // Recursively assign ancestor pointers
        // 
        void assign_ancestors(Node *node, Node **ancestors)
        {

            node->ancestors = new Node *[2*dim];
            for (size_t i = 0; i < 2*dim; ++i) {
                node->ancestors[i] = ancestors[i];
            }

            if (node->nodes) {

                Node **new_ancestors = new Node *[2*dim];

                for (size_t n = 0; n < nnodes; ++n) { 
                    if (node->nodes[n]) { 
                        for (size_t d = 0; d < dim; ++d) { 
                            if (n & (1 << d)) { 
                                new_ancestors[2*d] = node;
                                new_ancestors[2*d+1] = ancestors[2*d+1];
                            } else { 
                                new_ancestors[2*d] = ancestors[2*d];
                                new_ancestors[2*d+1] = node; 
                            }
                        }

                        assign_ancestors(node->nodes[n], new_ancestors);
                    }
                }

                delete[] new_ancestors;
            } 

        }


        void build_skiphierarchy(Node *level)
        { 
            //recursively build skip hierarchy until we get an empty level
            //should be O(log(n)) levels w.h.p.
            while (level) {
     
                level = build_skiplevel(level); 

                //store level and continue
                if (level) levels.push_back(level);
            } 
        }

        bool verify_parents(Node *node)
        {

            if (node->nodes) { 
                for (size_t n = 0; n < nnodes; ++n) { 
                    if (node->nodes[n]) {
                        if (node->nodes[n]->parent != node) return false;
                        verify_parents(node->nodes[n]); 
                    }
                }
            }

            return true;
        }

        bool verify_node_in_level(Node *root, Node *node)
        {
            if (root == node) return true;

            bool result = false;
            if (root->nodes) { 
                for (size_t n = 0; n < nnodes; ++n) { 
                    if (root->nodes[n]) {
                        result |= verify_node_in_level(root->nodes[n], node); 
                    }
                }
            }

            return result;

        }

        Node *build_skiplevel(Node *node)
        { 
            Node *copy = 0; 

            //hit a leaf -- nothing to do
            if (!node->nodes) return 0;

            //search through children for leaf nodes
            //we work at the branch node level so as to avoid creating nodes without children
            //and then deleting them
            for (size_t n = 0; n < nnodes; ++n) { 
                if (node->nodes[n]) {

                    //hit a leaf node
                    if (node->nodes[n]->nodes == 0) {
                        //decide with p=0.5 whether or not to keep the node
                        if ((double)rand() / (double)RAND_MAX < 0.5) { 

                            //if so, and we haven't made a copy yet, do so
                            if (!copy) copy = make_branchnode(node);

                            //create copy of existing point at proper location
                            copy->nodes[n] = new Node;
                            copy->nodes[n]->parent = copy;

                            //link nodes for skip hierarchy
                            node->nodes[n]->above = copy->nodes[n]; 
                            copy->nodes[n]->below = node->nodes[n];

                            //copy point in
                            copy->nodes[n]->pt = node->nodes[n]->pt;
                            for (size_t d = 0; d < dim; ++d) {
                                copy->nodes[n]->mid[d] = node->nodes[n]->mid[d];
                                copy->nodes[n]->radius = node->nodes[n]->radius; 
                            } 
                        }
                    } else {
                        //build level recursively
                        Node *child = build_skiplevel(node->nodes[n]);

                        //found a valid child, make copy if necessary and set pointer
                        if (child) {
                            if (!copy) copy = make_branchnode(node);
                            copy->nodes[n] = child;
                            copy->nodes[n]->parent = copy;

                            //link nodes for skip hierarchy
                            node->nodes[n]->above = copy->nodes[n]; 
                            copy->nodes[n]->below = node->nodes[n];
                        } 
                    }
                }
            }

            return copy;
        } 

        //make a new branch node and set pointers for skip hierarchy 
        Node *make_branchnode(Node *node)
        { 
            //make new branch node and zero out pointers
            Node *copy = new Node;
            copy->nodes = new Node *[nnodes];
            for (size_t n = 0; n < nnodes; ++n) {
                copy->nodes[n] = 0;
            } 

            //copy mid, side 
            for (size_t d = 0; d < dim; ++d) {
                copy->mid[d] = node->mid[d];
                copy->radius = node->radius; 
            }

            //link nodes for skip hierarchy
            node->above = copy;
            copy->below = node;

            return copy; 
        } 

        //
        // Search for smallest square q inside p that is equidistant to pt
        //
        Node *skip_search_equidistant(const Point &pt, Node *p, double p_dist)
        {

            //find highest point in skip hierarchy containing p, call it q
            Node *q_in_q0 = p;
            Node *last_q_in_q0 = 0;
            Node *q_in_qi = p;

            while (q_in_qi->above != 0) {
                q_in_qi = q_in_qi->above;
            }

            while (true) {

                // if we hit a leaf, we are done
                if (!q_in_q0->nodes) {
                    break;
                }

                //first check to see if q in Q0 has two equidistant children to p
                size_t nequidistant = 0;

                if (last_q_in_q0 != q_in_q0) {

                    for (size_t n = 0; n < nnodes; ++n) { 
                        //calculate distance to each of the non-zero children
                        if (q_in_q0->nodes[n]) {

                            double min_dist = min_pt_dist_to_node(pt, q_in_q0->nodes[n]);

                            //if close enough, count as equidistant
                            if (fabs(min_dist - p_dist) < EQUIDISTANT_EPSILON) {
                                ++nequidistant; 
                                if (nequidistant == 2) {
                                    break; 
                                }
                            }
                        } 
                    }

                    last_q_in_q0 = q_in_q0;
                }

                if (nequidistant == 2) break;

                //if not, see if q in Qi has an equidistant child to p
                size_t equidistant_child = -1;
                if (q_in_qi->nodes) {
                    for (size_t n = 0; n < nnodes; ++n) { 
                        if (q_in_qi->nodes[n]) {

                            double min_dist = min_pt_dist_to_node(pt, q_in_qi->nodes[n]);

                            //if close enough, count as equidistant
                            if (fabs(min_dist - p_dist) < EQUIDISTANT_EPSILON) {
                                equidistant_child = n;
                                break; 
                            }
                        }
                    }
                }

                //if equidistant child exists, move to it -- otherwise drop to lower level
                if (equidistant_child != -1) { 
                    q_in_qi = q_in_qi->nodes[equidistant_child];
                    q_in_q0 = q_in_q0->nodes[equidistant_child]; 
                } else {
                    if (q_in_qi->below) q_in_qi = q_in_qi->below;
                    else break;
                }
            }

            return q_in_q0;
        }

        double min_pt_dist_to_node(const Point &pt, Node *node)
        {
            bool inside = true; 
            double min_dist = std::numeric_limits<double>::max();
            for (size_t d = 0; d < dim; ++d) { 
        
                double dist; 
                if (pt[d] < node->mid[d] - node->radius) { 
                    dist = node->mid[d] - node->radius - pt[d];
                    inside = false;
                } else if (pt[d] > node->mid[d] + node->radius) {
                    dist = pt[d] - (node->mid[d] + node->radius); 
                    inside = false;
                }

                if (dist < min_dist) min_dist = dist; 
            } 

            if (inside) return 0.0;
            else return min_dist*min_dist;
        }

        double max_pt_dist_to_node(const Point &pt, Node *node)
        { 
            double max_dist = std::numeric_limits<double>::min();
            for (size_t d = 0; d < dim; ++d) { 
                double dist;
                dist = fabs(node->mid[d] - node->radius - pt[d]);
                if (dist > max_dist) max_dist = dist; 

                dist = fabs(pt[d] - node->mid[d] + node->radius);
                if (dist > max_dist) max_dist = dist; 
            } 

            return max_dist*max_dist;
        } 

        bool in_node(const Node *node, const Point &pt)
        {
            bool in = true;

            for (size_t d = 0; d < dim; ++d) { 
                if (root->mid[d] - root->radius - pt[d] > locate_eps || pt[d] - root->mid[d] - root->radius > locate_eps) { 
                    in = false; 
                    break; 
                } 
            } 

            return in; 
        } 

        void delete_worker(Node *node)
        {
            if (node->nodes) {
                for (size_t n = 0; n < nnodes; ++n) { 
                    if (node->nodes[n]) delete_worker(node->nodes[n]); 
                }

                delete node->nodes;
            }

            delete node->ancestors;
            delete node; 
        }    

};

#endif

