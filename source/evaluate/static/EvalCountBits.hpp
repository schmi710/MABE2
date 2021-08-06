/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalCountBits.hpp
 *  @brief MABE Evaluation module for counting the number of ones (or zeros) in an output.
 */

#ifndef MABE_EVAL_COUNT_BITS_H
#define MABE_EVAL_COUNT_BITS_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../core/Genome.hpp"


#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalCountBits : public Module {
  private:
    Collection target_collect;

    std::string bits_trait;
    std::string fitness_trait;
    bool count_type;   // =0 for counts zeros, or =1 for count ones.

  public:
    EvalCountBits(mabe::MABE & control,
                  const std::string & name="EvalCountBits",
                  const std::string & desc="Evaluate bitstrings by counting ones (or zeros).",
                  const std::string & _btrait="bits",
                  const std::string & _ftrait="fitness",
                  bool _ctype=1)
      : Module(control, name, desc)
      , target_collect(control.GetPopulation(0))
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , count_type(_ctype)
    {
      SetEvaluateMod(true);
    }
    ~EvalCountBits() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(count_type, "count_type", "Which type of bit should we count? (0 or 1)");
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

        // Count the number of ones in the bit sequence.
        //const emp::BitVector & bits = org.GetVar<emp::BitVector>(bits_trait);

        Genome::Head head = org.GetVar<Genome::Head>(bits_trait);

        double fitness = 0.0;
        double total = 0.0;

        while (head.IsValid()) {
          fitness += head.ReadBit();
          total += 1;
        }
        
        //(double) bits.CountOnes();

        // If we were supposed to count zeros, subtract ones count from total number of bits.
        if (count_type == 0) fitness = total - fitness;

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

  MABE_REGISTER_MODULE(EvalCountBits, "Evaluate bitstrings by counting ones (or zeros).");
}

#endif
