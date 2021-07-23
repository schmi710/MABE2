/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2020.
 *
 *  @file  BytesOrgGenomeHead.hpp
 *  @brief An organism consisting of a series of bytes.  Uses new Genome class and the head.
 *  @note Status: ALPHA
 */

#ifndef MABE_BYTES_GENOME_HEAD_ORGANISM_H
#define MABE_BYTES_GENOME_HEAD_ORGANISM_H

#include "../core/MABE.hpp"
#include "../core/Organism.hpp"
#include "../core/OrganismManager.hpp"
#include "../core/Genome.hpp"

#include "emp/bits/BitVector.hpp"
#include "emp/math/Distribution.hpp"
#include "emp/math/random_utils.hpp"

namespace mabe {

  class BytesOrgGenomeHead : public OrganismTemplate<BytesOrgGenomeHead> {
  protected:
    TypedGenome<std::byte> genome;

  public:
    BytesOrgGenomeHead(OrganismManager<BytesOrgGenomeHead> & _manager)
      : OrganismTemplate<BytesOrgGenomeHead>(_manager), genome(100) { }
    BytesOrgGenomeHead(const BytesOrgGenomeHead &) = default;
    BytesOrgGenomeHead(BytesOrgGenomeHead &&) = default;
    BytesOrgGenomeHead(const TypedGenome<std::byte> & in, OrganismManager<BytesOrgGenomeHead> & _manager)
      : OrganismTemplate<BytesOrgGenomeHead>(_manager), genome(in) { }
    BytesOrgGenomeHead(size_t N, OrganismManager<BytesOrgGenomeHead> & _manager)
      : OrganismTemplate<BytesOrgGenomeHead>(_manager), genome(N) { }
    ~BytesOrgGenomeHead() { ; }

    struct ManagerData : public Organism::ManagerData {
      double mut_prob = 0.01;            ///< Probability of each bit mutating on reproduction.
      std::string output_name = "bits";  ///< Name of trait that should be used to access bits.
      emp::Binomial mut_dist;            ///< Distribution of number of mutations to occur.
      emp::BitVector mut_sites;          ///< A pre-allocated vector for mutation sites. 
      bool init_random = true;           ///< Should we randomize ancestor?  (false = all zeros)
      double max_value = 1+0xFF;           ///< The maximum byte value (exclusive, not inclusive)
    };

    /// Use "to_string" to convert.
    std::string ToString() const override { return "using bytes: " + emp::to_string(genome); }

    size_t Mutate(emp::Random & random) override { 
      genome.SetTheDeets(SharedData().mut_prob, std::byte(0), SharedData().max_value);
      return genome.Mutate(random); 
      }

    void Randomize(emp::Random & random) override { genome.Randomize(random); }

    void Initialize(emp::Random & random) override { if (SharedData().init_random) genome.Randomize(random); }

    /// Put the bits in the correct output position.
    void GenerateOutput() override {
      SetVar<Genome::Head>(SharedData().output_name, Genome::Head(genome));
    }

    /// Setup this organism type to be able to load from config.
    void SetupConfig() override {
      GetManager().LinkFuns<size_t>([this](){ return genome.GetSize(); },
                       [this](const size_t & N){ return genome.Resize(N); },
                       "N", "Number of bytes in organism");
      GetManager().LinkVar(SharedData().mut_prob, "mut_prob",
                      "Probability of each byte mutating on reproduction.");
      GetManager().LinkVar(SharedData().max_value, "max_value",
                      "The maximum byte value.");
      GetManager().LinkVar(SharedData().output_name, "output_name",
                      "Name of variable to contain byte sequence.");
      GetManager().LinkVar(SharedData().init_random, "init_random",
                      "Should we randomize ancestor?  (0 = all zeros)");
    }

    /// Setup this organism type with the traits it need to track.
    void SetupModule() override {
      // Setup the mutation distribution.
      SharedData().mut_dist.Setup(SharedData().mut_prob, genome.GetSize());

      // Setup the default vector to indicate mutation positions.
      SharedData().mut_sites.Resize(genome.GetSize());

      // Setup the output trait.
      GetManager().AddSharedTrait(SharedData().output_name,
                                  "Byte output from organism.",
                                  Genome::Head(genome));
    }
  };

  MABE_REGISTER_ORG_TYPE(BytesOrgGenomeHead, "Organism consisting of a series of N bytes. Uses new Genome class and the head.");
}

#endif
