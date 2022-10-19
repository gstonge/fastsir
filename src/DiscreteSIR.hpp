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

#ifndef DISCRETESIR_HPP_
#define DISCRETESIR_HPP_

#include "BaseContagion.hpp"

namespace fastsir
{//start of namespace fastsir


//class to simulate SIS process on networks
class DiscreteSIR : public BaseContagion
{
public:
    //Constructor
    DiscreteSIR(const EdgeList& edge_list, double recovery_probability,
                const std::vector<double>& infection_probability);

    //Accessors
    double get_lifetime() const
        {return infected_node_set_.size() == 0 ?
            std::numeric_limits<double>::infinity() : 1.;}

    //Mutators
    void clear();

protected:
    //Members
    double recovery_probability_;
    std::vector<double> infection_probability_; //per node in group
    std::vector<double> infection_propensity_; //Poisson rate equiv
    sset::SamplableSet<Node> infection_event_set_;
    sset::SamplableSet<Node> recovery_event_set_;
    std::poisson_distribution<int> poisson_dist_;
    std::binomial_distribution<int> binomial_dist_;

    //utility functions
    inline double get_infection_propensity(Node node) const
        {return infection_propensity_.at(get_infected_degree(node));}

    inline void update_infection_propensity(Node node, const Event& event);

    inline void infect(Node node);
    inline void recover(Node node);
    inline std::vector<Event> next_step();
};

}//end of namespace fastsir

#endif /* DISCRETESIR_HPP_ */

