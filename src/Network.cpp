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

#include "Network.hpp"
#include <numeric>

using namespace std;

namespace fastsir
{//start of namespace fastsir

//Constructor of the class provided an edge list
Network::Network(const EdgeList& edge_list) :
	adjacency_list_(), nodes_(),
    min_degree_(0), max_degree_(0)
{
	size_t nb_nodes = 0;
	//Determine the number of nodes
	for (size_t i = 0; i < edge_list.size(); i++)
    {
    	if (edge_list[i].first > nb_nodes)
    	{
    		nb_nodes = edge_list[i].first;
    	}
    	if (edge_list[i].second > nb_nodes)
    	{
    		nb_nodes = edge_list[i].second;
    	}
    }
    nb_nodes += 1; //the label starts to 0 by convention

    //Initialize adjacency lists
    adjacency_list_ = AdjacencyList(nb_nodes, vector<Node>());
    nodes_ = vector<Node>(nb_nodes);
    iota(nodes_.begin(),nodes_.end(),0);

    for (auto & edge : edge_list)
    {
        adjacency_list_[edge.first].push_back(edge.second);
        adjacency_list_[edge.second].push_back(edge.first);
    }

    //Determine min and max degree
    for (Node node : nodes_)
    {
        if (node == 0)
        {
            min_degree_ = degree(node);
        }
        if (degree(node) < min_degree_)
        {
            min_degree_ = degree(node);
        }
        if (degree(node) > max_degree_)
        {
            max_degree_ = degree(node);
        }
    }
}

}//end of namespace fastsir
