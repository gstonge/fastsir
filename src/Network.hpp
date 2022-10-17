/*
 * MIT License
 *
 * Copyright (c) 2022 Guillaume St-Onge
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef NETWORK_HPP_
#define NETWORK_HPP_

#include <utility>
#include <vector>
#include <memory>

namespace fastsir
{//start of namespace fastsir

typedef unsigned int Node;
typedef std::vector<std::pair<Node,Node>> EdgeList;
typedef std::vector<std::vector<Node> > AdjacencyList;


//Structure representing an undirected network
class Network
{
public:
    //Constructor
    Network(const EdgeList& edge_list);

    //Accessors
    std::size_t min_degree() const
        {return min_degree_;}
    std::size_t max_degree() const
        {return max_degree_;}

    std::size_t degree(Node node) const
    	{return adjacency_list_[node].size();}

    std::size_t size() const
        {return adjacency_list_.size();}
    std::size_t number_of_nodes() const
        {return adjacency_list_.size();}

    const std::vector<Node>& adjacent_nodes(Node node) const
    	{return adjacency_list_[node];}

    const std::vector<Node>& nodes() const
        {return nodes_;}

private:
    //Members
    AdjacencyList adjacency_list_;
    std::vector<Node> nodes_;
    std::size_t min_degree_;
    std::size_t max_degree_;

};

}//end of namespace fastsir

#endif /* NETWORK_HPP_ */
