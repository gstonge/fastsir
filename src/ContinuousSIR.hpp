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

#ifndef CONTINUOUSSIR_HPP_
#define CONTINUOUSSIR_HPP_

#include "BaseContagion.hpp"

namespace fastsir
{//start of namespace fastsir


//class to simulate SIR process on networks
class ContinuousSIR : public BaseContagion
{
public:
    //Constructor
    ContinuousSIR(const EdgeList& edge_list, double recovery_rate,
                  const std::vector<double>& infection_rate);

    //Accessors
    double get_lifetime() const
        {return infected_node_set_.size() == 0 ?
            std::numeric_limits<double>::infinity() : 1./(get_total_rate());}

    //Mutators
    void clear();

protected:
    //Members
    double recovery_rate_;
    std::vector<double> infection_rate_;
    sset::SamplableSet<Node> infection_event_set_;
    sset::SamplableSet<Node> recovery_event_set_;

    //utility functions
    inline double get_infection_rate(Node node) const
        {return infection_rate_.at(get_infected_degree(node));}
    inline double get_total_rate() const
        {return infection_event_set_.total_weight() + recovery_rate_*recovery_event_set_.total_weight();}

    inline void update_infection_rate(Node node, const Event& event);

    inline void infect(Node node);
    inline void recover(Node node);
    inline std::vector<Event> next_step();
};

}//end of namespace fastsir

#endif /* CONTINUOUSSIR_HPP_ */

