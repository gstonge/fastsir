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

#ifndef BASECONTAGION_HPP_
#define BASECONTAGION_HPP_

#include "SamplableSet/SamplableSet.hpp"
#include <unordered_set>
#include <unordered_map>
#include "Network.hpp"
#include <iostream>

namespace fastsir
{//start of namespace fastsir
enum NodeState {S, I, R, COUNT};
const unsigned int STATECOUNT = static_cast<unsigned int>(NodeState::COUNT);
enum Action {RECOVERY,INFECTION};
//enum Actor {GROUP,NODE};

typedef double Time;
typedef unsigned int Infector;
typedef unsigned int Infectee;
typedef unsigned int Generation; //of the infector
typedef std::size_t InfectedDegree;
//typedef std::tuple<Actor,Action,Label> Event;
typedef std::pair<Node,Action> Event;
typedef std::tuple<Time,std::size_t,std::size_t,std::size_t> MacroState;
typedef std::tuple<Time,Generation,Infector,Infectee,InfectedDegree> Transmission;
typedef std::unordered_map<Node,std::size_t> InfectedNeighborPosition;


//abstract class with more functionality to avoid overlapp between classes
class BaseContagion
{
public:
    //Constructor
    BaseContagion(const EdgeList& edge_list);

    //Accessors
    std::size_t size() const
        {return network_.size();}
    const std::vector<NodeState>& get_state_vector() const
        {return state_vector_;}
    const std::unordered_set<Node>& get_infected_node_set() const
        {return infected_node_set_;}
    const Network& get_network() const
        {return network_;}
    double get_current_time() const
        {return current_time_;}
    std::size_t get_number_of_infected_nodes() const
        {return infected_node_set_.size();}
    InfectedDegree get_infected_degree(Node node) const
        {return infected_neighbors_vector_.at(node).size();}
    std::vector<MacroState> get_macro_state_vector() const
        {return macro_state_vector_;}
    std::vector<Transmission> get_transmission_tree() const
        {return transmission_tree_;}
    MacroState get_current_macro_state() const
        {return std::make_tuple(current_time_,
                                network_.size()-infected_node_set_.size()-recovered_node_set_.size(),
                                infected_node_set_.size(),
                                recovered_node_set_.size());}

    //Mutators
    void seed(unsigned int seed)
        {gen_.seed(seed);}
    void infect_fraction(double fraction);
    void infect_node_set(const std::unordered_set<Node>& node_set);

    void clear();
    void reset();

    void evolve(double period, bool save_transmission_tree, bool save_macro_state);


protected:
    //Members
    Network network_;
    std::vector<NodeState> state_vector_;
    std::unordered_set<Node> infected_node_set_;
    std::unordered_set<Node> recovered_node_set_;
    std::unordered_map<Node,Generation> infection_generation_;
    std::vector<std::vector<Node>> infected_neighbors_vector_;
    std::vector<InfectedNeighborPosition> infected_neighbor_position_vector_;
    std::vector<MacroState> macro_state_vector_;
    std::vector<Transmission> transmission_tree_;

    double current_time_;
    double last_event_time_;
    sset::RNGType &gen_;
    mutable std::uniform_real_distribution<double> random_01_;

    //utility functions
    Node random_infected_neighbor(Node node) const;
    inline void store_current_macro_state();
    inline void update_transmission_tree(const std::vector<Event>& event_vector);
    inline void apply_events(const std::vector<Event>& event_vector);

    //dummy functions because abstract virtual function breaks binding
    virtual double get_lifetime() const
        {return 1.;}
    virtual void infect(Node node) {};
    virtual void recover(Node node) {};
    virtual std::vector<Event> next_step()
        {return std::vector<Event>();}

};

}//end of namespace fastsir

#endif /* BASECONTAGION_HPP_ */
