/*
 * MIT License
 *
 * Copyright (c) 2021 Guillaume St-Onge
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

#include "BaseContagion.hpp"
#include <cmath>
#include <optional>
#include <utility>
#include <iostream>
#include <exception>

using namespace std;

namespace fastsir
{//start of namespace fastsir

//constructor of the class
BaseContagion::BaseContagion(const EdgeList& edge_list):
    network_(edge_list),
    state_vector_(network_.size(), S),
    infected_node_set_(),
    recovered_node_set_(),
    infection_generation_(),
    infected_neighbors_vector_(network_.size(), vector<Node>()),
    infected_neighbor_position_vector_(network_.size(), InfectedNeighborPosition()),
    current_time_(0),
    last_event_time_(0),
    gen_(sset::BaseSamplableSet::gen_),
    random_01_()
{
}

//get a random node of the particular state in the group
Node BaseContagion::random_infected_neighbor(Node node) const
{
    unsigned int index = floor(random_01_(gen_)*infected_neighbors_vector_.at(node).size());
    return (infected_neighbors_vector_.at(node)).at(index);
}

//store the current macro state
inline void BaseContagion::store_current_macro_state()
{

    macro_state_vector_.emplace_back(current_time_,network_.size()-infected_node_set_.size()-recovered_node_set_.size(),
            infected_node_set_.size(),recovered_node_set_.size());
}

//store the current macro state
inline void BaseContagion::update_transmission_tree(const std::vector<Event>& event_vector)
{
    for (const Event& event : event_vector)
    {
        Action action = event.second;
        if (action == INFECTION)
        {
            Node infectee = event.first;
            //pick a random infector and get generation
            Node infector = random_infected_neighbor(infectee);
            Generation generation = infection_generation_.at(infector);
            infection_generation_[infectee] = generation + 1;
            //get number of infected neighbors
            InfectedDegree infected_degree = get_infected_degree(infectee);

            transmission_tree_.emplace_back(current_time_,generation,infector,infectee,infected_degree);
        }
    }
}

inline void BaseContagion::apply_events(const vector<Event>& event_vector)
{
    for (const Event& event : event_vector)
    {
        Node node = event.first;
        Action action = event.second;
        if (action == INFECTION)
        {
            infect(node);
        }
        else if (action == RECOVERY)
        {
            recover(node);
        }
    }
}


//infect a fraction of the nodes initially
void BaseContagion::infect_fraction(double fraction)
{
    unsigned int number_of_infection = floor(network_.size()*fraction);
    Node node;
    unsigned int count = 0;
    while (count < number_of_infection)
    {
        node = floor(random_01_(gen_)*network_.size());
        if (state_vector_[node] == S)
        {
            infect(node);
            infection_generation_[node] = 0;
            count += 1;
        }
    }
}

//infect a certain set of of nodes
void BaseContagion::infect_node_set(const std::unordered_set<Node>& node_set)
{
    for (Node node : node_set)
    {
        if (state_vector_[node] == S)
        {
            infect(node);
            infection_generation_[node] = 0;
        }
    }
}

//clear the state; as if all node became susceptible at this time
void BaseContagion::clear()
{
    //recover nodes
    unordered_set<Node> infected_node_set_copy = infected_node_set_;
    for (Node node : infected_node_set_copy)
    {
        recover(node);
    }
    //put recovered nodes back to susceptible
    for (Node node : recovered_node_set_)
    {
        state_vector_[node] = S;
    }
    recovered_node_set_.clear();
    //other clear
    infection_generation_.clear();
}

//clear and reset the process to initial state at time 0 (and clear history)
//clear all measures as well
void BaseContagion::reset()
{
    clear();
    //reset transmission tree and macro state vector
    macro_state_vector_.clear();
    transmission_tree_.clear();
    current_time_ = 0;
    last_event_time_ = 0;
}

//perform the evolution of the process over a period of time and perform
//measures after each decorrelation time if needed
void BaseContagion::evolve(double period, bool save_transmission_tree, bool save_macro_state)
{
    if (save_macro_state and (macro_state_vector_.size() == 0))
    {
        store_current_macro_state();
    }
    double initial_time = current_time_;

    vector<Event> event_vector;
    while((last_event_time_ + get_lifetime() - initial_time <= period) and isfinite(get_lifetime()))
    {
        event_vector = next_step();
        //save transmission tree
        if (save_transmission_tree and event_vector.size() > 0)
        {
            update_transmission_tree(event_vector);
        }
        apply_events(event_vector);
        //store macro state
        if (save_macro_state)
        {
            store_current_macro_state();
        }
    }
    current_time_ = initial_time + period;
}

}//end of namespace fastsir
