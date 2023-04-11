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

#include "AnnealedContinuousSIR.hpp"
#include <optional>
#include <utility>
#include <iostream>
#include <exception>
#include <cmath>
#include <limits>

using namespace std;

namespace fastsir
{//start of namespace fastsir

//dummy edge list of a star graph
EdgeList dummy_edge_list(int number_of_nodes)
{
    EdgeList v;
    for (int i = 0; i < number_of_nodes; i++)
    {
        v.push_back(make_pair(0,i));
    }
    return v;
}

//to calculate binomial coefficient
int BinomialCoefficient(const int n, const int k) {
  vector<int> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

//constructor of the class
AnnealedContinuousSIR::AnnealedContinuousSIR(int number_of_nodes, int degree, double recovery_rate,
        const vector<double>& infection_rate):
    BaseContagion(dummy_edge_list(number_of_nodes)),
    degree_(degree),
    binomial_coeff_(degree+1,0),
    recovery_rate_(recovery_rate),
    infection_rate_(infection_rate),
    infection_event_set_(1.,1.),
    recovery_event_set_(1.,1.)
{
    //by default, we populate the infection event set with all susceptible nodes
    for (Node node = 0; node < size(); node++)
    {
        infection_event_set_.insert(node,1.);
    }
    //precalculate the binomial coefficients
    for (int i = 0; i <= degree; i++)
    {
        binomial_coeff_[i] = BinomialCoefficient(degree_, i);
    }
}

//update the infection rate of a neighbor node
inline void AnnealedContinuousSIR::update_meanfield_rate()
{
    double prevalence = double(get_number_of_infected_nodes())/size();
    double sum = 0.;
    for (int i = 0; i <= degree_ ; i++)
    {
        sum += binomial_coeff_[i]*pow(prevalence,i)*pow(1-prevalence,degree_-i)*get_infection_rate(i);
    }
    meanfield_infection_rate_ = sum;
}

//infect a node
inline void AnnealedContinuousSIR::infect(Node node)
{
    if (state_vector_[node] == S)
    {
        infection_event_set_.erase(node);
        state_vector_[node] = I;
        infected_node_set_.insert(node);

        //create a recovery event for the node
        recovery_event_set_.insert(node, 1.);

        //update the meanfield
        update_meanfield_rate();
    }
    else
    {
        throw runtime_error("Infection attempt: the node is not susceptible");
    }
}

//recover a node
inline void AnnealedContinuousSIR::recover(Node node)
{
    if (state_vector_[node] == I)
    {
        state_vector_[node] = R;
        infected_node_set_.erase(node);
        recovered_node_set_.insert(node);
        //erase the recovery event for the node
        recovery_event_set_.erase(node);

        //update the meanfield
        update_meanfield_rate();
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}


//advance the process to the next step by performing infection/recovery
inline vector<Event> AnnealedContinuousSIR::next_step()
{
    current_time_ = last_event_time_ + get_lifetime();
    last_event_time_ = current_time_;
    vector<Event> event_vector;

    //determine if the next event is infection or recovery
    if ((meanfield_infection_rate_*infection_event_set_.total_weight()/get_total_rate()) > random_01_(gen_))
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
void AnnealedContinuousSIR::clear()
{
    BaseContagion::clear();
    infection_event_set_.clear(); //to avoid numerical error accumulation
    recovery_event_set_.clear(); //to avoid numerical error accumulation
    //by default, we populate the infection event set with all susceptible nodes
    for (Node node = 0; node < size(); node++)
    {
        infection_event_set_.insert(node,1.);
    }
}



}//end of namespace fastsir
