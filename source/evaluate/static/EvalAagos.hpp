/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalCountBits.hpp
 *  @brief MABE Evaluation module for counting the number of ones (or zeros) in an output.
 */

#ifndef MABE_EVAL_AAGOS_H
#define MABE_EVAL_AAGOS_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"

#include "emp/datastructs/reference_vector.hpp"
#include "emp/base/vector.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/math/random_utils.hpp"

namespace mabe {

  class EvalAagos : public Module {
  private:
    Collection target_collect;

    std::string bits_trait;
    std::string fitness_trait;
    emp::vector<emp::BitVector> gene_targets;
    size_t gene_size;
    size_t change_frequency;
    size_t change_magnitude;
    emp::Random random;

  public:
    EvalAagos(mabe::MABE & control,
                  const std::string & name="EvalAagos",
                  const std::string & desc="Evaluate a genome by its proportion of bits matching the gene targets.")
      : Module(control, name, desc)
      , target_collect(control.GetPopulation(0))
      , bits_trait("bits")
      , fitness_trait("fitness")
      , gene_targets(8)
      , gene_size(8)
      , change_frequency(4)
      , change_magnitude(1)
      , random(control.GetRandom())

    {
      SetEvaluateMod(true);
    }
    ~EvalAagos() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
    }

    void SetupModule() override {
      AddRequiredTrait<emp::vector<emp::BitVector>>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "Proportion of matching bits fitness value", 0.0);

      // randomly initialize the gene targets
      for (size_t i = 0; i < gene_targets.size(); i++) {
          emp::BitVector gene(gene_size);
          emp::RandomizeBitVector(gene, random, 0.5);
          gene_targets[i] = gene;
      }
    }

    void OnUpdate(size_t update) override {
      emp_assert(control.GetNumPopulations() >= 1);

      // every 'change_frequency', mutate the environment
      if (update % change_frequency == 0) {
        mutateEnvironment();
      }

      // Loop through the population and evaluate each organism.
      double max_fitness = 0.0;
      emp::Ptr<Organism> max_org = nullptr;
      mabe::Collection alive_collect( target_collect.GetAlive() );
      for (Organism & org : alive_collect) {        
        // Make sure this organism has its bit sequence ready for us to access.
        org.GenerateOutput();
        // get genome from organism
        const emp::vector<emp::BitVector> & gene_set = org.GetVar<emp::vector<emp::BitVector>>(bits_trait);
        // start fitness at 0
        double fitness = 0.0;
        for (size_t i = 0; i < gene_targets.size(); i++) {
            double matching_bits = 0.0;
            for (size_t j = 0; j < gene_size; j++) {
                if (gene_set[i][j] == gene_targets[i][j]) {
                    matching_bits++;
                }
            }
            double gene_fitness = matching_bits / gene_size;
            fitness += gene_fitness;
        }
        fitness = (fitness / gene_targets.size()) * 100;

        // Store the count on the organism in the fitness trait.
        org.SetVar<double>(fitness_trait, fitness);

        if (fitness > max_fitness) {
          max_fitness = fitness;
        }
      }
      std::cout << "Max " << fitness_trait << " = " << max_fitness << std::endl;
    }

    // mutate 'change_magnitude' random bits from any gene in gene targets
    void mutateEnvironment() {
      emp::vector<emp::BitVector> change_sites(gene_targets.size()); // create an empty change site vector of bitvectors
      for (size_t i = 0; i < gene_targets.size(); i++) { // set each bitvector to all 0's
        emp::BitVector gene_change_site(gene_size);
        gene_change_site.Clear();
        change_sites[i] = gene_change_site;
      }

      // mutations in gene_targets may appear on the same gene multiple times, but the same bit will not be mutated more than once
      for(size_t i = 0; i < change_magnitude; i++) {
        const size_t gene_index = random.GetUInt(gene_targets.size()); // randomly pick a gene to mutate
        const size_t bit_index = random.GetUInt(gene_size); // randomly pick a bit in the gene to mutate
        if (change_sites[gene_index][bit_index]) { --i; continue; }  // duplicate position; try again.
        change_sites[i].Set(bit_index); // mark site as mutated
        gene_targets[gene_index].Toggle(bit_index); // mutate single bit
      }
    }
  };

  MABE_REGISTER_MODULE(EvalAagos, "Evaluate a genome by its proportion of bits matching the gene targets.");
}

#endif