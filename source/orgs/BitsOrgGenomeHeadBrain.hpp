/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2020.
 *
 *  @file  BitsOrgGenomeHeadBrain.hpp
 *  @brief An organism consisting of a series of bits.  Uses new Genome class and the head and the Brain.
 *  @note Status: ALPHA
 */

#ifndef MABE_BITS_GENOME_HEAD_BRAIN_ORGANISM_H
#define MABE_BITS_GENOME_HEAD_BRAIN_ORGANISM_H

#include "../core/MABE.hpp"
#include "../core/Organism.hpp"
#include "../core/OrganismManager.hpp"
#include "../core/Genome.hpp"
#include "../core/Brain.hpp"

#include "emp/bits/BitVector.hpp"
#include "emp/math/Distribution.hpp"
#include "emp/math/random_utils.hpp"

namespace mabe {

  class BitsOrgGenomeHeadBrain : public OrganismTemplate<BitsOrgGenomeHeadBrain> {
  protected:
    TypedGenome<bool> genome;
    BasicBrain<Genome::Head> brain;
    emp::DataMap map;
  public:
    BitsOrgGenomeHeadBrain(OrganismManager<BitsOrgGenomeHeadBrain> & _manager)
      : OrganismTemplate<BitsOrgGenomeHeadBrain>(_manager), genome(100), brain() { }
    BitsOrgGenomeHeadBrain(const BitsOrgGenomeHeadBrain &) = default;
    BitsOrgGenomeHeadBrain(BitsOrgGenomeHeadBrain &&) = default;
    BitsOrgGenomeHeadBrain(const TypedGenome<bool> & in, OrganismManager<BitsOrgGenomeHeadBrain> & _manager)
      : OrganismTemplate<BitsOrgGenomeHeadBrain>(_manager), genome(in), brain() { }
    BitsOrgGenomeHeadBrain(size_t N, OrganismManager<BitsOrgGenomeHeadBrain> & _manager)
      : OrganismTemplate<BitsOrgGenomeHeadBrain>(_manager), genome(N), brain() { }
    ~BitsOrgGenomeHeadBrain() { ; }

    struct ManagerData : public Organism::ManagerData {
      double mut_prob = 0.01;            ///< Probability of each bit mutating on reproduction.
      std::string output_name = "bits";  ///< Name of trait that should be used to access bits.
      emp::Binomial mut_dist;            ///< Distribution of number of mutations to occur.
      emp::BitVector mut_sites;          ///< A pre-allocated vector for mutation sites. 
      bool init_random = true;           ///< Should we randomize ancestor?  (false = all zeros)
    };

    /// Use "to_string" to convert.
    std::string ToString() const override { return "ima head out: " + emp::to_string(genome); }

    size_t Mutate(emp::Random & random) override { return genome.Mutate(random, SharedData().mut_prob); }

    void Randomize(emp::Random & random) override { genome.Randomize(random); }

    void Initialize(emp::Random & random) override { if (SharedData().init_random) genome.Randomize(random); }

    /// Put the bits in the correct output position.
    void GenerateOutput() override {
      map.Set<emp::Ptr<Organism>>("organism", emp::Ptr<Organism>(this));
      map.Set<Genome::Head>("output", Genome::Head(genome));
      brain.process(map);
    }

    /// Setup this organism type to be able to load from config.
    void SetupConfig() override {
      GetManager().LinkFuns<size_t>([this](){ return genome.GetSize(); },
                       [this](const size_t & N){ return genome.Resize(N); },
                       "N", "Number of bits in organism");
      GetManager().LinkVar(SharedData().mut_prob, "mut_prob",
                      "Probability of each bit mutating on reproduction.");
      GetManager().LinkVar(SharedData().output_name, "output_name",
                      "Name of variable to contain bit sequence.");
      GetManager().LinkVar(SharedData().init_random, "init_random",
                      "Should we randomize ancestor?  (0 = all zeros)");

      map.AddVar<Genome::Head>("output", Genome::Head(genome));
      map.AddVar<std::string>("output_name", SharedData().output_name);
      map.AddVar<emp::Ptr<Organism>>("organism", emp::Ptr<Organism>(this));
    }

    /// Setup this organism type with the traits it need to track.
    void SetupModule() override {
      // Setup the mutation distribution.
      SharedData().mut_dist.Setup(SharedData().mut_prob, genome.GetSize());

      // Setup the default vector to indicate mutation positions.
      SharedData().mut_sites.Resize(genome.GetSize());

      // Setup the output trait.
      GetManager().AddSharedTrait(SharedData().output_name,
                                  "Bitset output from organism.",
      //                            emp::BitVector(0));
                                  Genome::Head(genome));
    }
  };

  MABE_REGISTER_ORG_TYPE(BitsOrgGenomeHeadBrain, "Organism consisting of a series of N bits. Uses new Genome class and the head and the brain.");
}

#endif
