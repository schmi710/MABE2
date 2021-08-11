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
#include "emp/bits/BitVector.hpp"

namespace mabe {

  class Brain {
    protected:

    public:
        Brain() {};
        virtual ~Brain() { ; }
        Brain(const Brain &) = default;
        Brain & operator=(const Brain &) = default;
        virtual void process(emp::DataMap org_map) = 0;
  };

    // remember to put this in a different file later
  template <typename T> class BasicBrain : public Brain { // returns whatever's in output
    public:
        BasicBrain() {};
        BasicBrain(const BasicBrain &) = default;
        void process(emp::DataMap org_map) override {
            org_map.Get<emp::Ptr<Organism>>("organism")->SetVar<T>(org_map.Get<std::string>("output_name"), org_map.Get<T>("output"));
        }
        emp::Ptr<Brain> Clone() {
            return emp::NewPtr<BasicBrain<T>>(*this); 
        }
  };

// remember to put this in a different file later
  class DynamicBrain : public Brain { // all DynamicOrg brains should be derived from DynamicBrain
    protected:
        template <typename T>
        void SetVar(emp::Ptr<emp::DataMap> data_map, const std::string & name, const T & value) {
            if (data_map->HasName(name) == false) data_map->AddVar<T>(name, value);
            else data_map->Set<T>(name, value);
        }

        template <typename T>
        void SetVar(emp::Ptr<emp::DataMap> data_map, size_t id, const T & value) {
            emp_assert(data_map->HasID(id), id);
            data_map->Set<T>(id, value);
        }

        std::string genomes_string = "genomes";
        std::string output_string = "output";
        std::string output_name_string = "output_name";


        emp::vector<size_t> genome_indices;

    public:
        DynamicBrain(): genome_indices() {};
        DynamicBrain(emp::vector<size_t> gis): genome_indices(gis) {};
        DynamicBrain(const DynamicBrain &) = default;

        virtual void rebuild(emp::DataMap org_map) = 0;

        virtual void initializeGenomes(emp::DataMap org_map) = 0;

        virtual void process(emp::DataMap org_map) = 0;

        virtual emp::Ptr<DynamicBrain> Clone() = 0;

        emp::vector<emp::Ptr<Genome>> GetGenomes(emp::DataMap org_map) {
            return org_map.Get<emp::vector<emp::Ptr<Genome>>>("genomes");
        }

        emp::Ptr<emp::DataMap> GetOutput(emp::DataMap org_map) {
            return org_map.Get<emp::Ptr<emp::DataMap>>("output");
        }

        std::string GetOutputName(emp::DataMap org_map) {
            return org_map.Get<std::string>("output_name");
        }

        
  };



  // remember to put this in a different file later
  class BasicBrainHead : public DynamicBrain { // returns a head of the first genome

    public:
        BasicBrainHead(emp::vector<size_t> gis): DynamicBrain(gis) {};
        BasicBrainHead(const BasicBrainHead &) = default;
        void process(emp::DataMap org_map) override {

            SetVar<Genome::Head>(GetOutput(org_map), GetOutputName(org_map), Genome::Head(*(GetGenomes(org_map)[genome_indices[0]])));

        }
        emp::Ptr<DynamicBrain> Clone() override {
            return emp::NewPtr<BasicBrainHead>(*this); 
        }

        void rebuild(emp::DataMap org_map) override {};

        void initializeGenomes(emp::DataMap org_map) override {};
  };

    // remember to put this in a different file later
   class MarkovBrain : public DynamicBrain {
    protected:
        size_t num_inputs = 8; // as a whole
        size_t num_outputs = 2; // as a whole
        size_t num_gate_inputs = 3; // to its gates
        size_t num_gate_outputs = 3; // to its gates
        size_t num_memory = 5;


        //emp::BitVector memory;


    public:

        struct Gate {
            emp::vector<int> input_ids;      ///< Vector of the input ids (larger is memory)
            emp::vector<int> output_ids;     ///< Vector of the output ids (larger is memory)
            emp::BitVector table;            ///< Size is a power of 2, indexed into by size_t naturally

            Gate(emp::vector<int> in, emp::vector<int> out, emp::BitVector t): input_ids(in), output_ids(out), table(t) {};

            Gate() {};
            ~Gate() { ; }
            Gate(const Gate &) = default;
            Gate & operator=(const Gate &) = default;


            friend std::string to_string(Gate const& self) {
                return emp::ToString(self.input_ids) + emp::ToString(self.output_ids) + self.table.ToString();
            }

        };

        emp::vector<Gate> gates;

        //MarkovBrain(): memory(num_memory), gates() {};
        //MarkovBrain(size_t mem): num_memory(mem), memory(mem), gates() {};
        //MarkovBrain(size_t inputs, size_t outputs, size_t mem): num_inputs(inputs), num_outputs(outputs), num_memory(mem), memory(mem), gates() {};

        MarkovBrain(emp::vector<size_t> gis): DynamicBrain(gis), gates() {};
        MarkovBrain(size_t mem): num_memory(mem), gates() {};
        MarkovBrain(size_t inputs, size_t outputs, size_t mem): num_inputs(inputs), num_outputs(outputs), num_memory(mem), gates() {};

        MarkovBrain(const MarkovBrain &) = default;
        void rebuild(emp::DataMap org_map) override {

            Genome::Head head(*(GetGenomes(org_map)[genome_indices[0]]));


            head.Reset();
            emp::vector<size_t> start_codons;
            while (head.IsValid()) {
                size_t pos = head.MultiFindInt(emp::vector<int>({42, 213}));
                if (pos >= 0) {
                    start_codons.push_back(pos);
                }  
            }
            gates.clear();
            for (size_t i: start_codons) {
                head.SetPosition(i);
                head.Advance(2); // get past the start codon
                emp::vector<int> input_ids = head.MultiReadScaledInt(num_gate_inputs, 0, num_inputs + num_memory);
                emp::vector<int> output_ids = head.MultiReadScaledInt(num_gate_outputs, 0, num_inputs + num_memory);
                emp::BitVector table = head.MultiReadScaledBit(pow(2,num_gate_inputs) * num_gate_outputs);
            
                gates.push_back(Gate(input_ids, output_ids, table));
            }
        }


        void process(emp::DataMap org_map) override {
            /*
            emp::BitVector output_and_new_memory(num_outputs + num_memory);
            emp::BitVector input_and_old_memory(num_inputs + num_memory);
            input_and_old_memory.Import(org_map.Get<emp::BitVector>("input"), 0); 
            input_and_old_memory.Import(memory, num_inputs); 

            for (Gate gate: gates) {
                size_t pos = 0;
                for (size_t i = 0; i < num_gate_inputs; i ++) { // little endian
                    if (input_and_old_memory.Get(gate.table.Get(i))) {
                        pos += pow(2, i);
                    }
                }
                emp::BitVector result = gate.table.Export(num_gate_outputs, pos * num_gate_outputs);
                for (size_t i = 0; i < num_gate_outputs; i ++) {
                    if (result.Get(i)) {
                        output_and_new_memory.Set(i, 1);
                    }
                } 
            }

            org_map.Get<emp::Ptr<Organism>>("organism")->SetVar<emp::BitVector>(org_map.Get<std::string>("output_name"), output_and_new_memory.Export(num_outputs, 0));
            memory = output_and_new_memory.Export(num_memory, num_outputs);
            */

           

           SetVar<emp::vector<Gate>>(GetOutput(org_map), GetOutputName(org_map), gates);
        }

        emp::Ptr<DynamicBrain> Clone() override {
            return emp::NewPtr<MarkovBrain>(*this); 
        }

        void initializeGenomes(emp::DataMap org_map) override {}; // should have something here
  };

};


  

  

#endif
