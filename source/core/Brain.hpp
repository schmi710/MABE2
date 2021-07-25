/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file  Brain.hpp
 *  @brief Base brain representation for organisms.
 */

#ifndef MABE_BRAIN_HPP
#define MABE_BRAIN_HPP


#include "emp/data/DataMap.hpp"
#include "Organism.hpp"
#include "Genome.hpp"
#include "orgs/BitsOrgGenomeHeadBrain.hpp"

namespace mabe {

  class Brain {
    public:
        Brain() {};
        ~Brain() { ; }
        Brain(const Brain &) = default;
        Brain & operator=(const Brain &) = default;
        virtual void process(emp::DataMap org_map) = 0;
  };

  template <typename T> class BasicBrain : public Brain {
    public:
        BasicBrain() {};
        BasicBrain(const BasicBrain &) = default;
        void process(emp::DataMap org_map) override {
            org_map.Get<emp::Ptr<Organism>>("organism")->SetVar<T>(org_map.Get<std::string>("output_name"), org_map.Get<T>("output"));
        }
  };

};


  

  

#endif
