/*
 * MIT License
 *
 * Copyright (c) 2023 Guillaume St-Onge
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

#include "ContinuousSIR.hpp"
#include <optional>
#include <utility>
#include <iostream>
#include <exception>
#include <cmath>
#include <limits>

using namespace std;

namespace fastsir
{//start of namespace fastsir

//constructor of the class
ContinuousSIR::ContinuousSIR(const EdgeList& edge_list, double recovery_rate,
        const std::vector<double>& infection_rate):
    BaseContagion(edge_list),
    recovery_rate_(recovery_rate),
    infection_rate_(infection_rate),
    infection_event_set_(1.,1.),
    recovery_event_set_(1.,1.)
{
    //calculate Poisson rate equivalent for each rate
    double min = std::numeric_limits<double>::infinity();
    double max = 0;
    for (double rate : infection_rate)
    {
        if (rate > max)
        {
            max = rate;
        }
        if (rate < min and rate > 0)
        {
            min = rate;
        }
    }
    infection_event_set_ = sset::SamplableSet<Node>(min,max); //set true bounds
}

//update the infection rate of a neighbor node
inline void ContinuousSIR::update_infection_rate(Node node, const Event& event)
{
    Node other_node = event.first;
    Action action = event.second;

    vector<Node>& infected_neighbors = infected_neighbors_vector_[node];
    unordered_map<Node,size_t>& infected_neighbor_position =  infected_neighbor_position_vector_[node];

    size_t size_before = infected_neighbors.size();
    if (action == RECOVERY)
    {
        size_t position = infected_neighbor_position[other_node];

        swap(infected_neighbors[position],infected_neighbors.back());

        //also, update the position of the node in the back
        Node back_node = infected_neighbors[position];
        infected_neighbor_position[back_node] = position;
        //pop
        infected_neighbors.pop_back();
        infected_neighbor_position.erase(other_node);

    }
    else if (action == INFECTION)
    {
        infected_neighbor_position[other_node] = infected_neighbors.size();
        infected_neighbors.push_back(other_node);
    }

    //update event set with new rate
    double new_rate;
     new_rate = get_infection_rate(node);
    if (new_rate > 0)
    {
        infection_event_set_.set_weight(node,new_rate);
    }
    else
    {
        infection_event_set_.erase(node);
    }
}

//infect a node
inline void ContinuousSIR::infect(Node node)
{
    if (state_vector_[node] == S)
    {
        infection_event_set_.erase(node);
        state_vector_[node] = I;
        infected_node_set_.insert(node);
        Event event = make_pair(node,INFECTION);
        for (Node neighbor : network_.adjacent_nodes(node))
        {
            if (state_vector_[neighbor] == S)
            {
                update_infection_rate(neighbor, event);
            }
        }
        //create a recovery event for the node
        recovery_event_set_.insert(node, 1.);
        //clear infected neighbors
        infected_neighbors_vector_[node].clear();
        infected_neighbor_position_vector_[node].clear();
    }
    else
    {
        throw runtime_error("Infection attempt: the node is not susceptible");
    }
}

//recover a node
inline void ContinuousSIR::recover(Node node)
{
    if (state_vector_[node] == I)
    {
        state_vector_[node] = R;
        infected_node_set_.erase(node);
        recovered_node_set_.insert(node);
        Event event = make_pair(node,RECOVERY);
        for (Node neighbor : network_.adjacent_nodes(node))
        {
            if (state_vector_[neighbor] == S)
            {
                update_infection_rate(neighbor, event);
            }
        }
        //erase the recovery event for the node
        recovery_event_set_.erase(node);
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}


//advance the process to the next step by performing infection/recovery
inline vector<Event> ContinuousSIR::next_step()
{
    current_time_ = last_event_time_ + get_lifetime();
    last_event_time_ = current_time_;
    vector<Event> event_vector;

    //determine if the next event is infection or recovery
    if ((infection_event_set_.total_weight()/get_total_rate()) > random_01_(gen_))
    {
        //infection event
        Node node = (infection_event_set_.sample()).value().first;
        event_vector.emplace_back(node,INFECTION);
    }
    else
    {
        //recovery event
        Node node = (recovery_event_set_.sample()).value().first;
        event_vector.emplace_back(node,RECOVERY);
    }

    return event_vector;
}


//clear the state; as if all node became susceptible at this time
//clear all measures as well
//overload BaseContagion
void ContinuousSIR::clear()
{
    BaseContagion::clear();
    infection_event_set_.clear(); //to avoid numerical error accumulation
    recovery_event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace fastsir
