/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file  EvalSuperPermsHead.hpp
 *  @brief MABE Evaluation module for finding the diversity of permutations in an output.
 */

#ifndef MABE_EVAL_SUPER_PERMS_HEAD_H
#define MABE_EVAL_SUPER_PERMS_HEAD_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../core/Genome.hpp"


#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalSuperPermsHead : public Module {
  private:
    Collection target_collect;

    std::string bits_trait;
    std::string fitness_trait;
    int count_length;   // =0 for counts zeros, or =1 for count ones.

  public:
    EvalSuperPermsHead(mabe::MABE & control,
                  const std::string & name="EvalSuperPermsHead",
                  const std::string & desc="Evaluate bitstrings by finding the diversity of permutations.",
                  const std::string & _btrait="bits",
                  const std::string & _ftrait="fitness",
                  int _clength=4)
      : Module(control, name, desc)
      , target_collect(control.GetPopulation(0))
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , count_length(_clength)
    {
      SetEvaluateMod(true);
    }
    ~EvalSuperPermsHead() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(count_length, "count_length", "How long are the patterns we're looking for?");
    }

    void SetupModule() override {
      AddRequiredTrait<Genome::Head>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "All-ones fitness value", 0.0);
    }

    void OnUpdate(size_t /* update */) override {
      emp_assert(control.GetNumPopulations() >= 1);

      // Loop through the population and evaluate each organism.
      double max_fitness = 0.0;
      emp::Ptr<Organism> max_org = nullptr;
      mabe::Collection alive_collect( target_collect.GetAlive() );
      for (Organism & org : alive_collect) {        
        // Make sure this organism has its bit sequence ready for us to access.
        org.GenerateOutput();

        Genome::Head & head = org.GetVar<Genome::Head>(bits_trait);

        size_t num_perms = pow(2,count_length);
        std::vector<int> counts(num_perms, 0);

/*
        while (head.IsValid()) {
          const emp::BitVector & temp_bits = head.MultiReadBit(count_length);
          head.Retreat();
          if (head.IsValid()) {
            head.Retreat(count_length-2);
            int ind = 0;
            for (int j = 0; j < count_length; j++) {
              ind += temp_bits[j] * pow(2,count_length - 1 - j);
            }
            counts[ind] ++;
          } else {
            break;
          }
          
        }*/
        const emp::BitVector & temp_bits = head.MultiReadBit(count_length);
        int ind = 0;
        for (int j = 0; j < count_length; j++) {
          ind += temp_bits[j] * pow(2,count_length - 1 - j);
        }
        counts[ind] ++;


        while (head.IsValid()) {
          ind = ((2 * ind) % ((int)pow(2,count_length))) + ((int)head.ReadBit());
          counts[ind] ++;
        }



        int min_set = *std::min_element(std::begin(counts), std::end(counts));
        double fitness = (num_perms * 2) * min_set;
        for (size_t i = 0; i < num_perms; i ++) {
            if (counts[i] > min_set) {
              fitness ++;
            }
        }



        // Store the count on the organism in the fitness trait.
        org.SetVar<double>(fitness_trait, fitness);

        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
        }
      }

      std::cout << "Max " << fitness_trait << " = " << max_fitness << std::endl;
    }
  };

  MABE_REGISTER_MODULE(EvalSuperPermsHead, "Evaluate bitstrings by finding the diversity of permutations.");
}

#endif
