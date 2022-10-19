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

#include "DiscreteSIR.hpp"
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
DiscreteSIR::DiscreteSIR(const EdgeList& edge_list, double recovery_probability,
        const std::vector<double>& infection_probability):
    BaseContagion(edge_list),
    recovery_probability_(recovery_probability),
    infection_probability_(infection_probability),
    infection_propensity_(),
    infection_event_set_(1.,1.),
    recovery_event_set_(1.,1.),
    poisson_dist_(1.),
    binomial_dist_(1,recovery_probability)
{
    //calculate Poisson rate equivalent for each probability
    double min = std::numeric_limits<double>::infinity();
    double max = 0;
    double propensity = 0;
    for (double prob : infection_probability)
    {
        propensity = -log(1.-prob);
        infection_propensity_.push_back(propensity);
        if (propensity > max)
        {
            max = propensity;
        }
        if (propensity < min and propensity > 0)
        {
            min = propensity;
        }
    }
    infection_event_set_ = sset::SamplableSet<Node>(min,max); //set true bounds
}

//update the infection propensity of a neighbor node
inline void DiscreteSIR::update_infection_propensity(Node node, const Event& event)
{
    Node other_node = event.first;
    Action action = event.second;

    vector<Node>& infected_neighbors = infected_neighbors_vector_[node];
    unordered_map<Node,size_t>& infected_neighbor_position =  infected_neighbor_position_vector_[node];
    if (action == RECOVERY)
    {
        //cout << "========start RECOVERY========" << endl;
        //put other node in the back position if not already
        //cout << "before CRIT" << endl;
        size_t position = infected_neighbor_position[other_node];

        //cout << "after CRIT" << endl;
        //cout << "position " << position << endl;
        //cout << "size " << infected_neighbors.size() << endl;
        swap(infected_neighbors[position],infected_neighbors.back());
        //cout << "after swap" << endl;
        //also, update the position of the node in the back
        Node back_node = infected_neighbors[position];
        infected_neighbor_position[back_node] = position;
        //pop
        infected_neighbors.pop_back();
        infected_neighbor_position.erase(other_node);
        //cout << "========end RECOVERY=========" << endl;
    }
    else if (action == INFECTION)
    {
        infected_neighbor_position[other_node] = infected_neighbors.size();
        infected_neighbors.push_back(other_node);
    }


    //update event set with new propensity
    double new_propensity = get_infection_propensity(node);
    if (new_propensity > 0)
    {
        infection_event_set_.set_weight(node,new_propensity);
    }
    else
    {
        infection_event_set_.erase(node);
    }
}

//infect a node
inline void DiscreteSIR::infect(Node node)
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
                update_infection_propensity(neighbor, event);
            }
        }
        //create a recovery event for the node
        recovery_event_set_.insert(node, 1.);
    }
    else
    {
        throw runtime_error("Infection attempt: the node is not susceptible");
    }
}

//recover a node
inline void DiscreteSIR::recover(Node node)
{
    if (state_vector_[node] == I)
    {
        state_vector_[node] = R;
        infected_node_set_.erase(node);
        recovered_node_set_.insert(node);
        Event event = make_pair(node,RECOVERY);
        //cout << "before update" << endl;
        for (Node neighbor : network_.adjacent_nodes(node))
        {
            if (state_vector_[neighbor] == S)
            {
                update_infection_propensity(neighbor, event);
            }
        }
        //cout << "after update" << endl;
        //erase the recovery event for the node
        recovery_event_set_.erase(node);
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}


//advance the process to the next step by performing infection/recovery
inline vector<Event> DiscreteSIR::next_step()
{
    current_time_ = last_event_time_ + get_lifetime();
    last_event_time_ = current_time_;

    pair<Node,double> node_weight_pair;
    //get the number of recoveries and assign them
    //cout << "before sampling new recovered"  << endl;
    binomial_dist_ = binomial_distribution<int>(
            recovery_event_set_.size(),recovery_probability_);
    int nb_rec = binomial_dist_(gen_);
    unordered_set<Node> new_susceptible; //use a set to discard repetition
    for (int i = 0; i < nb_rec; i++)
    {
        node_weight_pair = (recovery_event_set_.sample()).value();
        new_susceptible.insert(node_weight_pair.first);
    }
    //cout << "after sampling new recovered" << endl;

    //get the number of infections and assign them
    //cout << "before sampling new infected"  << endl;
    poisson_dist_ = poisson_distribution<int>(
            infection_event_set_.total_weight());
    int nb_inf = poisson_dist_(gen_);
    unordered_set<Node> new_infected;
    for (int i = 0; i < nb_inf; i++)
    {
        node_weight_pair = (infection_event_set_.sample()).value();
        new_infected.insert(node_weight_pair.first);
    }
    //cout << "after sampling new infected"  << endl;

    //return vector of events
    vector<Event> event_vector;
    for (Node node : new_susceptible)
    {
        event_vector.emplace_back(node,RECOVERY);
    }
    for (Node node : new_infected)
    {
        event_vector.emplace_back(node,INFECTION);
    }
    return event_vector;
}


//clear the state; as if all node became susceptible at this time
//clear all measures as well
//overload BaseContagion
void DiscreteSIR::clear()
{
    BaseContagion::clear();
    infection_event_set_.clear(); //to avoid numerical error accumulation
    recovery_event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace fastsir
