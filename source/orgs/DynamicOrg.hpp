/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2020.
 *
 *  @file  DynamicOrg.hpp
 *  @brief An organism consisting of a series of bits.  Uses new Genome class.
 *  @note Status: ALPHA
 */

#ifndef MABE_DYNAMIC_ORGANISM_H
#define MABE_DYNAMIC_ORGANISM_H

#include "../core/MABE.hpp"
#include "../core/Organism.hpp"
#include "../core/OrganismManager.hpp"
#include "../core/Genome.hpp"
#include "../brains/AagosBrain.hpp"

#include "emp/bits/BitVector.hpp"
#include "emp/math/Distribution.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/tools/string_utils.hpp"

namespace mabe {

  class DynamicOrg : public OrganismTemplate<DynamicOrg> {
  protected:
  public:


    emp::vector<double> mut_probs;
    emp::DataMap map;
    
    
    emp::vector<emp::Ptr<Genome>> genome_vec;
    emp::vector<emp::Ptr<DynamicBrain>> brain_vec;


    size_t size;

  
    DynamicOrg(OrganismManager<DynamicOrg> & _manager)
      : OrganismTemplate<DynamicOrg>(_manager), size(150){ 
          ClearAll();
      }
    //DynamicOrg(const DynamicOrg &) = default;
    DynamicOrg(const DynamicOrg & in): OrganismTemplate<DynamicOrg>(in), mut_probs(in.mut_probs), map(in.map), size(in.size){
    //DynamicOrg(const DynamicOrg & in): OrganismTemplate<DynamicOrg>(OrganismManager<DynamicOrg>((OrganismManager<DynamicOrg> &) in.GetManager())), mut_probs(in.mut_probs), map(in.map), genome(in.genome){
        ClearAll();
        for (size_t i = 0; i < in.genome_vec.size(); i++) {
            genome_vec.push_back(in.genome_vec[i]->Clone());
        }
        for (size_t i = 0; i < in.brain_vec.size(); i++) {
            brain_vec.push_back(in.brain_vec[i]->Clone());
        }
    }
    //DynamicOrg(DynamicOrg &&) = default;
    DynamicOrg(DynamicOrg && in): OrganismTemplate<DynamicOrg>(in.GetManager()), mut_probs(in.mut_probs), map(in.map), size(in.size){
        ClearAll();
        for (size_t i = 0; i < in.genome_vec.size(); i++) {
            genome_vec.push_back(in.genome_vec[i]->Clone());
        }
        for (size_t i = 0; i < in.brain_vec.size(); i++) {
            brain_vec.push_back(in.brain_vec[i]->Clone());
        }
    }
    DynamicOrg(const emp::vector<emp::Ptr<Genome>> genomes_in, const emp::vector<emp::Ptr<DynamicBrain>> brains_in, size_t N,
    OrganismManager<DynamicOrg> & _manager)
      : OrganismTemplate<DynamicOrg>(_manager), size(N){ 
            ClearAll();
            for (size_t i = 0; i < genomes_in.size(); i++) {
                genome_vec.push_back(genomes_in[i]->Clone());
            }
            for (size_t i = 0; i < brains_in.size(); i++) {
               brain_vec.push_back(brains_in[i]->Clone());
            }
      }
    DynamicOrg(size_t N, OrganismManager<DynamicOrg> & _manager)
      : OrganismTemplate<DynamicOrg>(_manager), size(N){ 
          ClearAll();
      }
    ~DynamicOrg() { 
        for (size_t i = 0; i < genome_vec.size(); i++) {
            genome_vec[i].Delete();
        }

        for (size_t i = 0; i < brain_vec.size(); i++) {
            brain_vec[i].Delete();
        }
    }

    void ClearAll() {
        for (size_t i = 0; i < genome_vec.size(); i++) {
            genome_vec[i].Delete();
        }

        for (size_t i = 0; i < brain_vec.size(); i++) {
            brain_vec[i].Delete();
        }
        genome_vec.clear();
        brain_vec.clear();
    }

    struct ManagerData : public Organism::ManagerData {
      std::string genomes = "bits_genome";  ///< string variable holds comma separated genome types
      std::string brains = "markov_brain";             ///< string variable holds comma separated brain types
      double mut_prob = 0.01;            ///< Probability of each bit mutating on reproduction.
      std::string output_name = "bits";  ///< Name of trait that should be used to access bits.
      emp::Binomial mut_dist;            ///< Distribution of number of mutations to occur.
      emp::BitVector mut_sites;          ///< A pre-allocated vector for mutation sites. 
      bool init_random = true;           ///< Should we randomize ancestor?  (false = all zeros)
    };

    // concatenates the ToString() of each genome separated by a new-line character and returns it
    // example of aagos ToString() would look like:
    // 1001011010101001101110010010100101
    // 16 4 0 9 12 5 8 13
    std::string ToString() const override {
      std::string out_string = "";
      
      for (size_t i = 0; i < genome_vec.size(); i++) {
        out_string += genome_vec[i]->ToString();
        out_string += "\n";
      }
      
      return out_string;
    }

    size_t Mutate(emp::Random & random) override { 
        size_t temp = 0;
        for (size_t i = 0; i < genome_vec.size(); i++) {
            temp += genome_vec[i]->Mutate(random, SharedData().mut_prob); 
        }
        return temp;
    }

    void Randomize(emp::Random & random) override { 
        for (size_t i = 0; i < genome_vec.size(); i++) {
            genome_vec[i]->Randomize(random);
        }
    }

    void Initialize(emp::Random & random) override {
      ClearAll();
      initGenomes(random);
      initBrains(random);
      if (SharedData().init_random) Randomize(random);
    }

    void initGenomes(emp::Random & random) {
      emp::vector<std::string> split_genomes;
      emp::slice(SharedData().genomes, split_genomes, ',');
      for(size_t i = 0; i < split_genomes.size(); i++) {
        if (split_genomes[i].compare("bit_genome") == 0) {
          genome_vec.push_back(emp::NewPtr<TypedGenome<bool>>(TypedGenome<bool>(size)));
        }

        if (split_genomes[i].compare("int_genome") == 0) {
          genome_vec.push_back(emp::NewPtr<TypedGenome<int>>(TypedGenome<int>(size)));
        }

        if (split_genomes[i].compare("byte_genome") == 0) {
          genome_vec.push_back(emp::NewPtr<TypedGenome<std::byte>>(TypedGenome<std::byte>(size)));
        }

        if (split_genomes[i].compare("val_genome") == 0) {
          genome_vec.push_back(emp::NewPtr<TypedGenome<double>>(TypedGenome<double>(size)));
        }
      }
      
    }

    void initBrains(emp::Random & random) {
      emp::vector<std::string> split_brains;
      emp::slice(SharedData().brains, split_brains, ',');
      map.Set<emp::vector<emp::Ptr<Genome>>>("genomes", genome_vec);
      for(size_t i = 0; i < split_brains.size(); i++) {
        if (split_brains[i].compare("aagos_brain") == 0) {
          //brain_vec.push_back(emp::NewPtr<AagosBrain>(genome_vec));             // Fix this with brain configuration. The brain should be told which genomes to use.
        } else if (split_brains[i].compare("markov_brain") == 0) {
          emp::Ptr<MarkovBrain> markov_brain = emp::NewPtr<MarkovBrain>(emp::vector<size_t>({0})); // Fix this with brain configuration. The brain should be told which genomes to use.
          brain_vec.push_back(markov_brain);
        } else if (split_brains[i].compare("basic_brain_head") == 0) {
          brain_vec.push_back(emp::NewPtr<BasicBrainHead>(emp::vector<size_t>({0})));               // Fix this with brain configuration. The brain should be told which genomes to use.
        }
        brain_vec[i]->initializeGenomes(map);
      }
      
    }

    size_t DoResize(size_t N) {
        size = N;
        for (size_t i = 0; i < genome_vec.size(); i++) {
            genome_vec[i]->Resize(N); 
        }
        return N;
    }



    /// Put the bits in the correct output position.
    void GenerateOutput() override {
      emp::DataMap temp;
  
      map.Set<emp::vector<emp::Ptr<Genome>>>("genomes", genome_vec);
      map.Set<emp::Ptr<emp::DataMap>>("data_map", & GetDataMap());
      map.Set<emp::Ptr<emp::DataMap>>("output", & temp);
      for (size_t i = 0; i < brain_vec.size(); i++) {
        brain_vec[i]->rebuild(map);
        brain_vec[i]->process(map);
      }

      SetVar<emp::DataMap>(SharedData().output_name, temp);
    }

    /// Setup this organism type to be able to load from config.
    void SetupConfig() override {
        GetManager().LinkVar(SharedData().genomes, "genomes", "Types of genomes to use.");
      GetManager().LinkVar(SharedData().brains, "brains", "Types of brains to use.");
      GetManager().LinkFuns<size_t>([this](){ return size; },
                       [this](const size_t & N){ return DoResize(N); },
                       "N", "Number of bits in organism");
      GetManager().LinkVar(SharedData().mut_prob, "mut_prob",
                      "Probability of each bit mutating on reproduction.");
      GetManager().LinkVar(SharedData().output_name, "output_name",
                      "Name of variable to contain bit sequence.");
      GetManager().LinkVar(SharedData().init_random, "init_random",
                      "Should we randomize ancestor?  (0 = all zeros)");


      map.AddVar<std::string>("output_name", SharedData().output_name);
      map.AddVar<emp::vector<emp::Ptr<Genome>>>("genomes", genome_vec);
      map.AddVar<emp::Ptr<emp::DataMap>>("data_map", & GetDataMap());
      emp::DataMap temp;
      map.AddVar<emp::Ptr<emp::DataMap>>("output", & temp);
    }

    /// Setup this organism type with the traits it need to track.
    void SetupModule() override {
      // Setup the mutation distribution.
      SharedData().mut_dist.Setup(SharedData().mut_prob, size);

      // Setup the default vector to indicate mutation positions.
      SharedData().mut_sites.Resize(size);

      // Setup the output trait.
      //TypedGenome<bool> temp(50);
      emp::DataMap temp;
      GetManager().AddSharedTrait(SharedData().output_name,
                                  "Bitset output from organism.",
                                  temp);
    }
  };

  MABE_REGISTER_ORG_TYPE(DynamicOrg, "A dynamically created organism composed of a series of Brains and Genomes.");
}

#endif