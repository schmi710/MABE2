/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file  AagosBrain.hpp
 *  @brief Aagos brain for use with AagosOrgGenome.hpp.
 */

#ifndef MABE_AAGOS_BRAIN_HPP
#define MABE_AAGOS_BRAIN_HPP

#include "../core/Brain.hpp"
#include "../core/Organism.hpp"
#include "emp/data/DataMap.hpp"
#include "../core/Genome.hpp"

#include "emp/base/vector.hpp"
#include "emp/bits/BitVector.hpp"

namespace mabe {
  class AagosBrain : public Brain {
    protected:
        emp::vector<emp::Ptr<Genome>> genomes;
    public:
      AagosBrain() {};
      AagosBrain(const AagosBrain &) = default;
      AagosBrain(const emp::vector<emp::Ptr<Genome>> & genomes_in) : genomes(genomes_in) { }
      void process(emp::DataMap org_map) override {
        size_t gene_size = 8;
        emp::vector<emp::BitVector> genes(genomes[0]->GetSize());
        for(size_t i = 0; i < genomes[1]->GetSize(); i++) {
          emp::BitVector gene(gene_size); // create empty bit vector which will hold separated genes
          gene.Clear(); // this function apparently only works for doubles right now, so I will manually set each bit to 0

          size_t current_position = genomes[0]->ReadInt(i); // start current position at genes start position
          for (size_t j = 0; j < gene_size; j++) {
            if (genomes[1]->ReadBit(current_position)) { // copy bit out of bits at current_position into gene at j
              gene.Toggle(j);
            }
            current_position++; // increase current position in genome
            if(current_position == genomes[1]->GetSize()) { // if current_position goes out of range, reset to 0
                current_position = 0;
            } 
          }
          genes[i] = gene; // insert constructed gene into output vector
        }
        org_map.Get<emp::Ptr<Organism>>("organism")->SetVar<emp::vector<emp::BitVector>>(org_map.Get<std::string>("output_name"), genes);
      }
  };
  //Karl Arne Johnson 
  //Marcelino Nosa Telmo
};

#endif