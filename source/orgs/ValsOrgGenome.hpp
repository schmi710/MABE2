/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file ValsOrgGenome.hpp
 *  @brief An organism consisting of a series of values of type double.  Uses Genome.
 *  @note Status: ALPHA
 */

#ifndef MABE_VALS_GENOME_ORGANISM_H
#define MABE_VALS_GENOME_ORGANISM_H

#include "../core/MABE.hpp"
#include "../core/Organism.hpp"
#include "../core/OrganismManager.hpp"
#include "../core/Genome.hpp"

#include "emp/base/vector.hpp"
#include "emp/math/Distribution.hpp"
#include "emp/math/random_utils.hpp"

namespace mabe {

  class ValsOrgGenome : public OrganismTemplate<ValsOrgGenome> {
  protected:
    TypedGenome<double> genome;
    double total = 0.0;        // Dynamic total of values in organism.

    // How do we enforce limits on values?
    enum BoundType {
      LIMIT_NONE=0,  // No boundary limit.  (e.g., in a 0 to 100 range, 103 would stay 103)
      LIMIT_CLAMP,   // Hard stop at boundary and stay there.    (e.g., 103 would go to 100)
      LIMIT_WRAP,    // Loop around through opposite boundary?   (e.g., 103 would go to 3)
      LIMIT_REBOUND, // Reflect back by amount limit was crossed (e.g., 103 would go to 97)
      LIMIT_ERROR    // Invalid limit type.
    };

    void CalculateTotal() {
      //for (double x : vals) total += x;
      for (double x : genome.GetAllData()) total += x;
      SetVar<double>(SharedData().total_name, total);
    }

  public:
    struct ManagerData : public Organism::ManagerData {
      std::string output_name = "vals";  ///< Name of trait that should be used to access values.
      std::string total_name = "total";  ///< Name of trait that indicate total of all values.
      double mut_prob = 0.01;            ///< Probability of position mutating on reproduction.
      double mut_size = 1.0;             ///< Standard deviation of mutations.
      double min_value = 0.0;            ///< Smallest that values are allowed to be.
      double max_value = 10000.0;          ///< Largest that values are allowed to be.
      BoundType upper_bound = LIMIT_REBOUND;
      BoundType lower_bound = LIMIT_REBOUND;

      // Helper member variables.
      emp::Binomial mut_dist;            ///< Distribution of number of mutations to occur.
      emp::BitVector mut_sites;          ///< A pre-allocated vector for mutation sites. 
      bool init_random = true;           ///< Should we randomize ancestor?  (false = all 0.0)

      // Helper functions.
      inline void ApplyBounds(double & value);              ///< Put a single value back in range.
      inline void ApplyBounds(emp::vector<double> & vals);  ///< Put all values back in range.
    };

    ValsOrgGenome(OrganismManager<ValsOrgGenome> & _manager)
      : OrganismTemplate<ValsOrgGenome>(_manager), genome(100), total(0.0) { }
    ValsOrgGenome(const ValsOrgGenome &) = default;
    ValsOrgGenome(ValsOrgGenome &&) = default;
    ValsOrgGenome(const TypedGenome<double> & in, OrganismManager<ValsOrgGenome> & _manager)
      : OrganismTemplate<ValsOrgGenome>(_manager), genome(in)
    {
      emp::vector<double> temp = genome.GetAllData();
      SharedData().ApplyBounds( temp );  // Make sure all data is within range.
      CalculateTotal();
    }
    ValsOrgGenome(size_t N, OrganismManager<ValsOrgGenome> & _manager)
      : OrganismTemplate<ValsOrgGenome>(_manager), genome(N), total(0.0) { }
    ~ValsOrgGenome() { ; }

    /// Use "to_string" to convert.
    std::string ToString() const override { return emp::to_string(genome, ":(TOTAL blah=", total, ")"); }

    size_t Mutate(emp::Random & random) override {
      // Identify number of and positions for mutations.
      /*const size_t num_muts = SharedData().mut_dist.PickRandom(random);
      emp::BitVector & mut_sites = SharedData().mut_sites;
      mut_sites.ChooseRandom(random, num_muts);

      // Trigger mutations at the identified positions.
      int mut_pos = mut_sites.FindOne();
      while (mut_pos != -1) {
        double & cur_val = vals[mut_pos];        // Identify the next site to mutate.
        total -= cur_val;                        // Remove old value from the total.
        cur_val += random.GetRandNormal();       // Mutate the value at the site.
        SharedData().ApplyBounds(cur_val);       // Make sure the value stays in the allowed range.
        total += cur_val;                        // Add the update value back into the total.
        mut_pos = mut_sites.FindOne(mut_pos+1);  // Move on to the next site to mutate.
      }*/
      genome.SetTheDeets(SharedData().mut_prob, SharedData().min_value, SharedData().max_value-SharedData().min_value);
      const size_t num_muts = genome.Mutate(random);

      total = 0.0;
      for (double x : genome.GetAllData()) total += x;

      SetVar<double>(SharedData().total_name, total);  // Store total in data map.
/*
      for (int i = 0; i < 20; i ++) {
        genome.WriteScaledByte(i, std::byte(i), std::byte(0), std::byte(40));
      }

      for (int i = 20; i < 30; i ++) {
        genome.WriteScaledByte(i, std::byte(i), std::byte(0), std::byte(30));
      }

      for (int i = 30; i < 40; i ++) {
        genome.WriteScaledBit(i, i % 2);
      }

      for (int i = 40; i < 50; i ++) {
        genome.WriteScaledByte(i, std::byte(i), std::byte(40), std::byte(50));
      }

      TypedGenome<bool> temp_genome(20);

      for (int i = 0; i < 20; i ++) {
        temp_genome.WriteScaledByte(i, std::byte(i), std::byte(0), std::byte(20));
      }

      for (int i = 50; i < 70; i ++) {
        genome.WriteScaledBit(i, temp_genome.ReadBit(i-50));
      }
*/

      return num_muts;
    }

    void Randomize(emp::Random & random) override {
      total = 0.0;
      //for (double & x : vals) {
      //  x = random.GetDouble(SharedData().min_value, SharedData().max_value);
      //  total += x;
      //}
      genome.SetTheDeets(SharedData().mut_prob, SharedData().min_value, SharedData().max_value-SharedData().min_value);
      genome.Randomize(random);
      for (double x : genome.GetAllData()) total += x;
      SetVar<double>(SharedData().total_name, total);  // Store total in data map.
    }

    void Initialize(emp::Random & random) override {
      if (SharedData().init_random) Randomize(random);
      //else { total = 0.0; for (double & x : vals) x = 0.0; }
      else { total = 0.0; genome.Clear(); }
    }


    /// Put the values in the correct output positions.
    void GenerateOutput() override {
      //SetVar<emp::vector<double>>(SharedData().output_name, vals);
      SetVar<emp::vector<double>>(SharedData().output_name, genome.GetAllData());
      SetVar<double>(SharedData().total_name, total);
    }

    /// Setup this organism type to be able to load from config.
    void SetupConfig() override {
      //GetManager().LinkFuns<size_t>([this](){ return vals.size(); },
      //                 [this](const size_t & N){ return vals.resize(N, 0.0); },
      //                 "N", "Number of values in organism");
      GetManager().LinkFuns<size_t>([this](){ return genome.GetSize(); },
                       [this](const size_t & N){ return genome.Resize(N, 0.0); },
                       "N", "Number of values in organism");
      GetManager().LinkVar(SharedData().mut_prob, "mut_prob",
                      "Probability of each value mutating on reproduction.");
      GetManager().LinkVar(SharedData().mut_size, "mut_size",
                      "Standard deviation on size of mutations.");
      GetManager().LinkVar(SharedData().min_value, "min_value",
                      "Lower limit for value fields.");
      GetManager().LinkVar(SharedData().max_value, "max_value",
                      "Upper limit for value fields.");
      GetManager().LinkMenu(
        SharedData().lower_bound, "lower_bound", "How should the lower limit be enforced?",
        LIMIT_NONE, "no_limit", "Allow values to be arbirarily low.",
        LIMIT_CLAMP, "clamp", "Reduce too-low values to min_value.",
        LIMIT_WRAP, "wrap", "Make low values loop around to maximum.",
        LIMIT_REBOUND, "rebound", "Make low values 'bounce' back up." );
      GetManager().LinkMenu(
        SharedData().upper_bound, "upper_bound", "How should the upper limit be enforced?",
        LIMIT_NONE, "no_limit", "Allow values to be arbirarily high.",
        LIMIT_CLAMP, "clamp", "Reduce too-high values to max_value.",
        LIMIT_WRAP, "wrap", "Make high values loop around to minimum.",
        LIMIT_REBOUND, "rebound", "Make high values 'bounce' back down." );
      GetManager().LinkVar(SharedData().output_name, "output_name",
                      "Name of variable to contain set of values.");
      GetManager().LinkVar(SharedData().total_name, "total_name",
                      "Name of variable to contain total of all values.");
      GetManager().LinkVar(SharedData().init_random, "init_random",
                      "Should we randomize ancestor?  (0 = all 0.0)");
    }

    /// Setup this organism type with the traits it need to track.
    void SetupModule() override {
      // Setup the mutation distribution.
      //SharedData().mut_dist.Setup(SharedData().mut_prob, vals.size());
      SharedData().mut_dist.Setup(SharedData().mut_prob, genome.GetSize());

      // Setup the default vector to indicate mutation positions.
      //SharedData().mut_sites.Resize(vals.size());
      SharedData().mut_sites.Resize( genome.GetSize());

      // Setup the output trait.
      GetManager().AddSharedTrait(SharedData().output_name,
                                  "Value vector output from organism.",
                                  //emp::vector<double>(vals.size()));
                                  emp::vector<double>(genome.GetSize()));
      // Setup the output trait.
      GetManager().AddSharedTrait(SharedData().total_name,
                                  "Total of all organism outputs.",
                                  0.0);
    }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////
  //  Helper functions....

  void ValsOrgGenome::ManagerData::ApplyBounds(double & value) {
    if (value > max_value) {
      switch (upper_bound) {
        case LIMIT_NONE:    break;
        case LIMIT_CLAMP:   value = max_value; break;
        case LIMIT_WRAP:    value -= (max_value - min_value); break;
        case LIMIT_REBOUND: value = 2 * max_value - value; break;
        case LIMIT_ERROR:   break;  // For now; perhaps do somethign with error?
      }
    }
    else if (value < min_value) {
      switch (lower_bound) {
        case LIMIT_NONE:    break;
        case LIMIT_CLAMP:   value = min_value; break;
        case LIMIT_WRAP:    value += (max_value - min_value); break;
        case LIMIT_REBOUND: value = 2 * min_value - value; break;
        case LIMIT_ERROR:   break;  // For now; perhaps do somethign with error?
      }
    }
  }

  void ValsOrgGenome::ManagerData::ApplyBounds(emp::vector<double> & vals) {
    const size_t range_size = max_value - min_value;

    switch (upper_bound) {
    case LIMIT_NONE: break;
    case LIMIT_CLAMP:
      for (double & value : vals) {
        if (value > max_value) value = max_value;
        else if (value < min_value) value = min_value;        
      }
      break;
    case LIMIT_WRAP:
      for (double & value : vals) {
        if (value > max_value) value -= range_size;
        else if (value < min_value) value += range_size;        
      }
      break;
    case LIMIT_REBOUND:
      for (double & value : vals) {
        if (value > max_value) value = 2 * max_value - value;
        else if (value < min_value) value = 2 * min_value - value;
      }
      break;
    default:
      break; // Should probably throw error.
    }
  }


  MABE_REGISTER_ORG_TYPE(ValsOrgGenome, "Organism consisting of a series of N floating-point values. Uses Genome.");
}

#endif
