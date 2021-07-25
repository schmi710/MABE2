/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file  Genome.hpp
 *  @brief Base genome representation for organisms.
 */

#ifndef MABE_GENOME_HPP
#define MABE_GENOME_HPP

#include <limits>

#include "emp/base/error.hpp"
#include "emp/math/Random.hpp"
#include "emp/meta/TypeID.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/math/Distribution.hpp"
#include "emp/math/random_utils.hpp"

namespace mabe {

  // Interface class for all genome types.
  class Genome {
  public:
    // Note: No constructors given; Genomes are pure virtual and as such must always   
    //       be constructed via a derived class.

    virtual ~Genome() { ; }

    virtual emp::Ptr<Genome> Clone() = 0;         // Make an exact copy of this genome.
    virtual emp::Ptr<Genome> CloneProtocol() = 0; // Copy everything in this genome except sequence.

    virtual size_t GetSize() const = 0;           // Return the size of this genome using underlying type.
    virtual void Resize(size_t new_size) = 0;     // Set new size using underlying type.
    virtual void Resize(size_t new_size, double default_val) = 0;
    virtual size_t GetNumBytes() const = 0;       // Return number of bytes of data in this genome.
    virtual void SetSizeRange(size_t min_size, size_t max_size) = 0; // Put limits on genome size.

    virtual bool IsValid(size_t pos) const { return pos < GetSize(); }
    virtual size_t ValidatePosition(size_t pos) const { return pos; } // If circular, should mod, etc.

    virtual void Randomize(emp::Random &, size_t /*pos*/) = 0; // Randomize only at one locus

    // Randomize whole genome
    virtual void Randomize(emp::Random & random) {
      for (size_t i = 0; i < GetSize(); i++) Randomize(random, i); 
    }
    // Test for mutations in whole genome; return number of mutations occurred.
    virtual size_t Mutate(emp::Random & random) = 0;

    // Human-readable (if not easily understandable) shorthand representations.
    virtual std::string ToString() const { return "[unknown]"; }
    virtual void FromString(std::string & in) { emp_error("Cannot read genome from the string \"" + in + "\""); }

    // Potentially more compressed or structured formats for saving/loading genomes.
    // TODO:
    // Archive & Serialize(Archive & archive) = 0;

    // Genome accessors for individual values…
    virtual int ReadInt(size_t index) const = 0;
    virtual double ReadDouble(size_t index) const = 0;
    virtual std::byte ReadByte(size_t index) const = 0;
    virtual bool ReadBit(size_t index) const = 0;

    virtual emp::vector<int> MultiReadInt(size_t index, size_t dist) const = 0;
    virtual emp::vector<double> MultiReadDouble(size_t index, size_t dist) const = 0;
    virtual emp::vector<std::byte> MultiReadByte(size_t index, size_t dist) const = 0;
    virtual emp::BitVector MultiReadBit(size_t index, size_t dist) const = 0;

    virtual void WriteInt(size_t index, int value) =  0;
    virtual void WriteDouble(size_t index, double value) =  0;
    virtual void WriteByte(size_t index, std::byte value) =  0;
    virtual void WriteBit(size_t index, bool value) =  0;

    virtual void MultiWriteInt(size_t index, emp::vector<int> values) = 0;
    virtual void MultiWriteDouble(size_t index, emp::vector<double> values) = 0;
    virtual void MultiWriteByte(size_t index, emp::vector<std::byte> values) = 0;
    virtual void MultiWriteBit(size_t index, emp::BitVector values) = 0;


    virtual int ReadScaledInt(size_t index, int min_v, int max_v) const = 0;
    virtual double ReadScaledDouble(size_t index, double min_v, double max_v) const = 0;
    virtual std::byte ReadScaledByte(size_t index, std::byte min_v, std::byte max_v) const = 0;
    virtual bool ReadScaledBit(size_t index, bool min_v, bool max_v) const = 0;

    virtual emp::vector<int> MultiReadScaledInt(size_t index, size_t dist, int min_v, int max_v) const = 0;
    virtual emp::vector<double> MultiReadScaledDouble(size_t index, size_t dist, double min_v, double max_v) const = 0;
    virtual emp::vector<std::byte> MultiReadScaledByte(size_t index, size_t dist, std::byte min_v, std::byte max_v) const = 0;
    virtual emp::BitVector MultiReadScaledBit(size_t index, size_t dist, bool min_v, bool max_v) const = 0;

    virtual void WriteScaledInt(size_t index, int value, int min_v, int max_v) =  0;
    virtual void WriteScaledDouble(size_t index, double value, double min_v, double max_v) =  0;
    virtual void WriteScaledByte(size_t index, std::byte value, std::byte min_v, std::byte max_v) =  0;
    virtual void WriteScaledBit(size_t index, bool value, bool min_v, bool max_v) =  0;

    virtual void MultiWriteScaledInt(size_t index, emp::vector<int> values, int min_v, int max_v) = 0;
    virtual void MultiWriteScaledDouble(size_t index, emp::vector<double> values, double min_v, double max_v) = 0;
    virtual void MultiWriteScaledByte(size_t index, emp::vector<std::byte> values, std::byte min_v, std::byte max_v) = 0;
    virtual void MultiWriteScaledBit(size_t index, emp::BitVector values, bool min_v, bool max_v) = 0;

    virtual void Clear() = 0;
    virtual void SetAll() = 0;

    // NEED READS AND WRITES THAT ARE RESTRICTED TO A RANGE
    // EXAMPLE:
    // virtual int ReadInt(size_t index, int min, int max) const = 0;

    class Head {
    protected:
      emp::Ptr<Genome> genome;  // Attached genome.
      size_t pos;               // What position is this head located at?
      int direction = 1;        // Direction forward is forward.

      static const constexpr unsigned int NORMAL=0;
      static const constexpr unsigned int END_OF_GENOME=1;
      static const constexpr unsigned int END_OF_CHROMOSOME=2;
      static const constexpr unsigned int INVALID=4;

      unsigned int state = NORMAL;

    public:
      Head() {};
      Head(Genome & in_genome, size_t in_pos=0, int in_dir=1)
        : genome(&in_genome), pos(in_pos), direction(in_dir) { }
      Head(const Head &) = default;
      Head & operator=(const Head &) = default;

      // Check for Heads at the same position and same direction!
      bool operator==(const Head & in) const {
          return genome==in.genome && pos==in.pos && direction==in.direction;
      }
      bool operator!=(const Head & in) const { return !(*this == in); }

      bool IsValid() const { return genome->IsValid(pos); }

      // @CAO VALIDATIONS SHOULD HAPPEN BEFORE ACTIONS, NOT AFTER.

      Head & SetPosition(size_t in_pos) {
        pos = genome->ValidatePosition(in_pos);
        if (!IsValid()) state = INVALID;
        return *this;
      }

      Head & operator+ (const size_t dist) const
      {
        Head newHead(*this);
        return newHead.Advance(dist);
      }

      Head & operator- (const size_t dist) const
      {
        Head newHead(*this);
        return newHead.Retreat(dist);
      }

      Head & operator+= (const size_t dist)
      {
        return this->Advance(dist);
      }

      Head & operator-= (const size_t dist)
      {
        return this->Retreat(dist);
      }

      Head & Advance(const size_t factor=1) { return SetPosition(pos + direction * factor); }
      Head & Retreat(const size_t factor=1) { return Advance(-factor); }

      int ReadInt() { int out = IsValid() ? genome->ReadInt(pos) : 0; Advance(); return out; }
      double ReadDouble() { double out = IsValid() ? genome->ReadDouble(pos) : 0; Advance(); return out; }
      std::byte ReadByte() { std::byte out = IsValid() ? genome->ReadByte(pos) : std::byte(0); Advance(); return out; }
      bool ReadBit() { bool out = IsValid() ? genome->ReadBit(pos) : 0; Advance(); return out; }

      Head & WriteInt(int value) { if (IsValid()) genome->WriteInt(pos, value); return Advance(); }
      Head & WriteDouble(double value) { if (IsValid()) genome->WriteDouble(pos, value); return Advance(); }
      Head & WriteByte(std::byte value) { if (IsValid()) genome->WriteByte(pos, value); return Advance(); }
      Head & WriteBit(bool value) { if (IsValid()) genome->WriteBit(pos, value); return Advance(); }

      emp::vector<int> MultiReadInt(size_t dist) { 
        emp::vector<int> outvec = genome->MultiReadInt(pos, dist);
        Advance(dist); 
        return outvec;
      }
      emp::vector<double> MultiReadDouble(size_t dist) { 
        emp::vector<double> outvec = genome->MultiReadDouble(pos, dist);
        Advance(dist); 
        return outvec;
      }
      emp::vector<std::byte> MultiReadByte(size_t dist) { 
        emp::vector<std::byte> outvec = genome->MultiReadByte(pos, dist);
        Advance(dist); 
        return outvec;
      }
      emp::BitVector MultiReadBit(size_t dist) { 
        emp::BitVector outvec = genome->MultiReadBit(pos, dist);
        Advance(dist); 
        return outvec;
      }

      Head & MultiWriteInt(emp::vector<int> values) { 
        genome->MultiWriteInt(pos, values);
        return Advance(values.size()); 
      }
      Head & MultiWriteDouble(emp::vector<double> values) {
        genome->MultiWriteDouble(pos, values);
        return Advance(values.size()); 
      }
      Head & MultiWriteByte(emp::vector<std::byte> values) {
        genome->MultiWriteByte(pos, values);
        return Advance(values.size()); 
      }
      Head & MultiWriteBit(emp::BitVector values) {
        genome->MultiWriteBit(pos, values);
        return Advance(values.GetSize()); 
      }

      // The head will end after the target
      size_t FindInt(int target) {
        while (IsValid()) {
          if (target == genome->ReadInt(pos)) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t MultiFindInt(emp::vector<int> target) {
        while (IsValid()) {
          if (target == genome->MultiReadInt(pos, target.size())) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t FindIntRange(int target_low, int target_high) {
        while (IsValid()) {
          int val = genome->ReadInt(pos);
          if (target_low <= val && target_high >= val) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t MultiFindIntRange(emp::vector<int> target_low, emp::vector<int> target_high) {
        while (IsValid()) {
          emp::vector<int> val = genome->MultiReadInt(pos, target_low.size());
          bool check = true;
          for (size_t i = 0; i < target_low.size(); i ++) {
            if (target_low[i] > val[i] || target_high[i] < val[i]) {
              check = false;
              break;
            }
          }
          if (check) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t FindDouble(double target) {
        while (IsValid()) {
          if (target == genome->ReadDouble(pos)) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t MultiFindDouble(emp::vector<double> target) {
        while (IsValid()) {
          if (target == genome->MultiReadDouble(pos, target.size())) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t FindDoubleRange(double target_low, double target_high) {
        while (IsValid()) {
          double val = genome->ReadDouble(pos);
          if (target_low <= val && target_high >= val) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t MultiFindDoubleRange(emp::vector<double> target_low, emp::vector<double> target_high) {
        while (IsValid()) {
          emp::vector<double> val = genome->MultiReadDouble(pos, target_low.size());
          bool check = true;
          for (size_t i = 0; i < target_low.size(); i ++) {
            if (target_low[i] > val[i] || target_high[i] < val[i]) {
              check = false;
              break;
            }
          }
          if (check) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t FindByte(std::byte target) {
        while (IsValid()) {
          if (target == genome->ReadByte(pos)) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t MultiFindByte(emp::vector<std::byte> target) {
        while (IsValid()) {
          if (target == genome->MultiReadByte(pos, target.size())) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t FindByteRange(std::byte target_low, std::byte target_high) {
        while (IsValid()) {
          std::byte val = genome->ReadByte(pos);
          if (target_low <= val && target_high >= val) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t MultiFindByteRange(emp::vector<std::byte> target_low, emp::vector<std::byte> target_high) {
        while (IsValid()) {
          emp::vector<std::byte> val = genome->MultiReadByte(pos, target_low.size());
          bool check = true;
          for (size_t i = 0; i < target_low.size(); i ++) {
            if (target_low[i] > val[i] || target_high[i] < val[i]) {
              check = false;
              break;
            }
          }
          if (check) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t FindBit(bool target) {
        while (IsValid()) {
          if (target == genome->ReadBit(pos)) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target
      size_t MultiFindBit(emp::BitVector target) {
        while (IsValid()) {
          if (target == genome->MultiReadBit(pos, target.GetSize())) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t FindBitRange(bool target_low, bool target_high) {
        if (target_low) return FindBit(true);
        if (!target_high) return FindBit(false);
        Advance();
        return pos;
      }

      // The head will end after the target, note that the bounds are both inclusive
      size_t MultiFindBitRange(emp::BitVector target_low, emp::BitVector target_high) {
        while (IsValid()) {
          emp::BitVector val = genome->MultiReadBit(pos, target_low.size());
          bool check = true;
          for (size_t i = 0; i < target_low.GetSize(); i ++) {
            if (target_low[i] > val[i] || target_high[i] < val[i]) {
              check = false;
              break;
            }
          }
          if (check) {
            Advance();
            return pos;
          }
          Advance();
        }
        return -1;
      }

      // NEED READS AND WRITES THAT ARE RESTRICTED TO A RANGE
      // EXAMPLE OF RANGED-READ:
      //int ReadInt(int min, int max) {
      //  int out = IsValid() ? genome->ReadInt(pos, min, max) : 0; Advance();
      //  return out;
      //}

      void Reset() { pos = 0; state = NORMAL; direction = 1; }
      void ReverseDirection() { direction = -direction; }

      bool AtBegin() { return pos == 0; }
      bool AtEnd() { return pos == genome->GetSize(); }

      void Randomize(emp::Random & random) { genome->Randomize(random, pos); }

    };

    Head GetHead(size_t position=0, int dir=1) { return Head(*this, position, dir); }
    Head begin() { return GetHead(); }
    Head end() { return GetHead(GetSize()); }
    Head rbegin() { return GetHead(GetSize(), -1); }
    Head rend() { return GetHead(0, -1); }
  };


  template <typename LOCUS_T>
  class TypedGenome : public Genome {
  protected:
    using locus_t = LOCUS_T;
    using this_t = TypedGenome<LOCUS_T>;

    emp::vector<locus_t> data;                        // Actual data in the genome.
    double mut_p = 0.0;                               // Mutation probability (LOTS TO DO HERE!)
    size_t min_size = 0;
    size_t max_size = static_cast<size_t>(-1);

    double alphabet_size = 4.0;

    emp::Binomial mut_dist;
    emp::BitVector mut_sites;

    locus_t min_val = static_cast<locus_t>(0);

  public:
    TypedGenome<LOCUS_T>() { }
    TypedGenome<LOCUS_T>(size_t _size): data(_size) { 
      min_size = _size;
    }
    TypedGenome<LOCUS_T>(const this_t &) = default;


 



    emp::Ptr<Genome> Clone() override { return emp::NewPtr<this_t>(*this); }
    emp::Ptr<Genome> CloneProtocol() override {
      emp::Ptr<this_t> out_ptr = emp::NewPtr<this_t>();
      out_ptr->mut_p = mut_p;
      out_ptr->min_size = min_size;
      out_ptr->max_size = max_size;
      return out_ptr;
    }

    size_t GetSize() const override { return data.size(); }
    void Resize(size_t new_size) override { data.resize(new_size); }
    void Resize(size_t new_size, double default_val) override {
      data.resize(new_size, static_cast<locus_t>(default_val));
    }
    size_t GetNumBytes() const override { return sizeof(locus_t) * GetSize(); }
    void SetSizeRange(size_t _min, size_t _max) override { min_size = _min, max_size = _max; }

    void Randomize(emp::Random & random) override {
      for (size_t i = 0; i < GetSize(); i++) Randomize(random, i); 
    }

    void Randomize(emp::Random & random, size_t pos) override {
      data[pos] = static_cast<locus_t>( min_size + random.GetDouble(alphabet_size) );
    }

    // Human-readable (if not easily understandable) shorthand representations.
    // @CAO... needs to be done properly!
    std::string ToString() const override { 
      
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        std::string s = "";
        for (size_t i = 0; i < GetSize(); i++) s += std::to_string((int) data[i]) + " ";
        return s;
      } else {
        return emp::to_string(data); 
      }
      
    }
    void FromString(std::string & in) override { emp_error("Cannot read genome from the string \"" + in + "\""); }

    // Potentially more compressed or structured formats for saving/loading genomes.
    // TODO:
    // Archive & Serialize(Archive & archive) = 0;

    // Genome accessors for individual values…
    int ReadInt(size_t index) const override { return static_cast<int>(data[index]); }
    double ReadDouble(size_t index) const override { return static_cast<double>(data[index]); }
    std::byte ReadByte(size_t index) const override { return static_cast<std::byte>(data[index]); }
    bool ReadBit(size_t index) const override { return static_cast<bool>(data[index]); }


    emp::vector<int> MultiReadInt(size_t index, size_t dist) const override {
      emp::vector<int> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadInt(index + i) : 0; 
      }
      return outvec;
    }
    emp::vector<double> MultiReadDouble(size_t index, size_t dist) const override {
      emp::vector<double> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadDouble(index + i) : 0; 
      }
      return outvec;
    }
    emp::vector<std::byte> MultiReadByte(size_t index, size_t dist) const override {
      emp::vector<std::byte> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadByte(index + i) : std::byte(0); 
      }
      return outvec;
    }
    emp::BitVector MultiReadBit(size_t index, size_t dist) const override {
      emp::BitVector outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadBit(index + i) : 0; 
      }
      return outvec;
    }

    void WriteInt(size_t index, int value) override { data[index] = static_cast<locus_t>(value); }
    void WriteDouble(size_t index, double value) override { data[index] = static_cast<locus_t>(value); }
    void WriteByte(size_t index, std::byte value) override { data[index] = static_cast<locus_t>(value); }
    void WriteBit(size_t index, bool value) override { data[index] = static_cast<locus_t>(value); }

    void MultiWriteInt(size_t index, emp::vector<int> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteInt(index + i, values[i]); 
      }
    }
    void MultiWriteDouble(size_t index, emp::vector<double> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteDouble(index + i, values[i]); 
      }
    }
    void MultiWriteByte(size_t index, emp::vector<std::byte> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteByte(index + i, values[i]); 
      }
    }
    void MultiWriteBit(size_t index, emp::BitVector values) override {
      size_t dist = values.GetSize();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteBit(index + i, values[i]); 
      }
    }

    // Scaled reads and writes: converts value to double, subtracts the minimum value, divides by the range, then multiplies by alphabet size and adds minimum value
    int ReadScaledInt(size_t index, int min_v, int max_v) const override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        return min_v + static_cast<int>(static_cast<double>(static_cast<int>(data[index]) - static_cast<int>(min_val))/alphabet_size*(max_v-min_v)); 
      } else {
        return min_v + static_cast<int>(static_cast<double>(data[index] - min_val)/alphabet_size*(max_v-min_v)); 
      }
    }
    double ReadScaledDouble(size_t index, double min_v, double max_v) const override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        return min_v + static_cast<double>(static_cast<int>(data[index]) - static_cast<int>(min_val))/alphabet_size*(max_v-min_v); 
      } else {
        return min_v + static_cast<double>(data[index] - min_val)/alphabet_size*(max_v-min_v); 
      }
    }
    std::byte ReadScaledByte(size_t index, std::byte min_v, std::byte max_v) const override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        return static_cast<std::byte>(static_cast<int>(min_v) + static_cast<int>(static_cast<double>(static_cast<int>(data[index]) - static_cast<int>(min_val))/alphabet_size*(static_cast<int>(max_v)-static_cast<int>(min_v)))); 
      } else {
        return static_cast<std::byte>(static_cast<int>(min_v) + static_cast<int>(static_cast<double>(data[index] - min_val)/alphabet_size*(static_cast<int>(max_v)-static_cast<int>(min_v))));
      } 
    }
    bool ReadScaledBit(size_t index, bool min_v=0, bool max_v=1) const override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        return min_v + static_cast<bool>(static_cast<double>(static_cast<int>(data[index]) - static_cast<int>(min_val))/alphabet_size*(max_v-min_v));  
      } else {
        return min_v + static_cast<bool>(static_cast<double>(data[index] - min_val)/alphabet_size*(max_v-min_v)); 
      }
    }
    emp::vector<int> MultiReadScaledInt(size_t index, size_t dist, int min_v, int max_v) const override {
      emp::vector<int> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledInt(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::vector<double> MultiReadScaledDouble(size_t index, size_t dist, double min_v, double max_v) const override {
      emp::vector<double> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledDouble(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::vector<std::byte> MultiReadScaledByte(size_t index, size_t dist, std::byte min_v, std::byte max_v) const override {
      emp::vector<std::byte> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledByte(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::BitVector MultiReadScaledBit(size_t index, size_t dist, bool min_v=0, bool max_v=1) const override {
      emp::BitVector outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledBit(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    void WriteScaledInt(size_t index, int value, int min_v, int max_v) override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        data[index] = static_cast<locus_t>(static_cast<int>(min_val) + static_cast<int>(static_cast<double>(value-min_v)*alphabet_size/(max_v-min_v))); 
      } else {
        data[index] = min_val + static_cast<locus_t>(static_cast<double>(value-min_v)*alphabet_size/(max_v-min_v)); 
      }
    }
    void WriteScaledDouble(size_t index, double value, double min_v, double max_v) override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        data[index] = static_cast<locus_t>(static_cast<int>(min_val) + static_cast<int>((value-min_v)*alphabet_size/(max_v-min_v))); 
      } else {
        data[index] = min_val + static_cast<locus_t>((value-min_v)*alphabet_size/(max_v-min_v)); 
      }
    }
    void WriteScaledByte(size_t index, std::byte value, std::byte min_v, std::byte max_v) override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        data[index] = static_cast<locus_t>(static_cast<int>(min_val) + static_cast<int>(static_cast<double>(static_cast<int>(value)-static_cast<int>(min_v))*alphabet_size/(static_cast<int>(max_v)-static_cast<int>(min_v))));
      } else {
        data[index] = min_val + static_cast<locus_t>(static_cast<double>(static_cast<int>(value)-static_cast<int>(min_v))*alphabet_size/(static_cast<int>(max_v)-static_cast<int>(min_v)));
      }
    }
    void WriteScaledBit(size_t index, bool value, bool min_v=0, bool max_v=1) override { 
      if constexpr (sizeof(locus_t) <= sizeof(int)) {
        data[index] =static_cast<locus_t>( static_cast<int>(min_val) + static_cast<int>(static_cast<double>(value-min_v)*alphabet_size/(max_v-min_v))); 
      } else {
        data[index] = min_val + static_cast<locus_t>(static_cast<double>(value-min_v)*alphabet_size/(max_v-min_v)); 
      }
    }
    void MultiWriteScaledInt(size_t index, emp::vector<int> values, int min_v, int max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledInt(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledDouble(size_t index, emp::vector<double> values, double min_v, double max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledDouble(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledByte(size_t index, emp::vector<std::byte> values, std::byte min_v, std::byte max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledByte(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledBit(size_t index, emp::BitVector values, bool min_v=0, bool max_v=1) override {
      size_t dist = values.GetSize();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledBit(index + i, values[i], min_v, max_v); 
      }
    }
    

    size_t Mutate(emp::Random & random) override {
      
      mut_dist.Setup(mut_p, data.size()); // NEEDED HERE?????
      mut_sites.Resize(data.size()); // NEEDED HERE?????

      const size_t num_muts = mut_dist.PickRandom(random);

        
      if (num_muts == 0) return 0;
      if (num_muts == 1) {
        const size_t pos = random.GetUInt(data.size());
        WriteDouble(pos, static_cast<double>(min_val) + random.GetDouble(alphabet_size) );
        return 1;
      }

      // Only remaining option is num_muts > 1.
      mut_sites.Clear();
      for (size_t i = 0; i < num_muts; i++) {
        const size_t pos = random.GetUInt(data.size());
        if (mut_sites[pos]) { --i; continue; }  // Duplicate position; try again.
        mut_sites.Set(pos);
        WriteDouble(pos, static_cast<double>(min_val) + random.GetDouble(alphabet_size) );
      }

      return num_muts;
    }

    void Clear() override { 
      for (size_t i = 0; i < data.size(); i++) {
        WriteDouble(i, 0 );
      }
    }

    void SetAll() override { 
      for (size_t i = 0; i < data.size(); i++) {
        WriteDouble(i, static_cast<double>(min_val) + alphabet_size );
      }
    }

    emp::vector<locus_t> GetAllData() const {
      return data;
    }

    void SetTheDeets(double _mut, locus_t _min_val, double _alphabet_size) { // remove this once genomemanager exists
      mut_p = _mut;
      min_val = _min_val;
      alphabet_size = _alphabet_size;
    };

    
 
  };

  template <>
  class TypedGenome<bool> : public Genome {
    protected:
    using this_t = TypedGenome<bool>;

    emp::BitVector data;                        // Actual data in the genome.
    double mut_p = 0.0;                               // Mutation probability (LOTS TO DO HERE!)
    size_t min_size = 0;
    size_t max_size = static_cast<size_t>(-1);

    double alphabet_size = 2.0;

    emp::Binomial mut_dist;
    emp::BitVector mut_sites;

  public:
    TypedGenome<bool>() { }
    TypedGenome<bool>(size_t _size): data(_size) { 
      min_size = _size;
    }
    TypedGenome<bool>(const this_t &) = default;


 



    emp::Ptr<Genome> Clone() override { return emp::NewPtr<this_t>(*this); }
    emp::Ptr<Genome> CloneProtocol() override {
      emp::Ptr<this_t> out_ptr = emp::NewPtr<this_t>();
      out_ptr->mut_p = mut_p;
      out_ptr->min_size = min_size;
      out_ptr->max_size = max_size;
      return out_ptr;
    }

    size_t GetSize() const override { return data.size(); }
    void Resize(size_t new_size) override { data.resize(new_size); }
    void Resize(size_t new_size, double default_val) override {
      data.resize(new_size);
      if (static_cast<bool>(default_val)) {
        data.SetAll();
      } else {
        data.Clear();
      }
    }
    size_t GetNumBytes() const override { return sizeof(bool) * GetSize(); }
    void SetSizeRange(size_t _min, size_t _max) override { min_size = _min, max_size = _max; }

    void Randomize(emp::Random & random) override {
      emp::RandomizeBitVector(data, random, 0.5);
    }

    void Randomize(emp::Random & random, size_t pos) override {
      data[pos] = static_cast<bool>( random.GetDouble(alphabet_size) );
    }

    // Human-readable (if not easily understandable) shorthand representations.
    // @CAO... needs to be done properly!
    std::string ToString() const override { return data.ToString(); }
    void FromString(std::string & in) override { data = emp::BitVector(in); }

    // Potentially more compressed or structured formats for saving/loading genomes.
    // TODO:
    // Archive & Serialize(Archive & archive) = 0;

    // Genome accessors for individual values…
    int ReadInt(size_t index) const override { return static_cast<int>(data[index]); }
    double ReadDouble(size_t index) const override { return static_cast<double>(data[index]); }
    std::byte ReadByte(size_t index) const override { return static_cast<std::byte>(data[index]); }
    bool ReadBit(size_t index) const override { return static_cast<bool>(data[index]); }

    emp::vector<int> MultiReadInt(size_t index, size_t dist) const override {
      emp::vector<int> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadInt(index + i) : 0; 
      }
      return outvec;
    }
    emp::vector<double> MultiReadDouble(size_t index, size_t dist) const override{
      emp::vector<double> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadDouble(index + i) : 0; 
      }
      return outvec;
    }
    emp::vector<std::byte> MultiReadByte(size_t index, size_t dist) const override{
      emp::vector<std::byte> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadByte(index + i) : std::byte(0); 
      }
      return outvec;
    }
    emp::BitVector MultiReadBit(size_t index, size_t dist) const override{
      return data.Export(dist, index);
    }

    void WriteInt(size_t index, int value) override { data.Set(index,static_cast<bool>(value)); }
    void WriteDouble(size_t index, double value) override { data.Set(index, static_cast<bool>(value)); }
    void WriteByte(size_t index, std::byte value) override { data.Set(index, static_cast<bool>(value)); }
    void WriteBit(size_t index, bool value) override { data.Set(index, value); }

    void MultiWriteInt(size_t index, emp::vector<int> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteInt(index + i, values[i]); 
      }
    }
    void MultiWriteDouble(size_t index, emp::vector<double> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteDouble(index + i, values[i]); 
      }
    }
    void MultiWriteByte(size_t index, emp::vector<std::byte> values) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteByte(index + i, values[i]); 
      }
    }
    void MultiWriteBit(size_t index, emp::BitVector values) override {
      data = data.Import( values, index );
    }


    // Scaled reads and writes: converts value to double, subtracts the minimum value, divides by the range, then multiplies by alphabet size and adds minimum value
    int ReadScaledInt(size_t index, int min_v, int max_v) const override { return min_v + static_cast<int>(static_cast<double>(data[index])/alphabet_size*(max_v-min_v)); }
    double ReadScaledDouble(size_t index, double min_v, double max_v) const override { return min_v + static_cast<double>(data[index])/alphabet_size*(max_v-min_v); }
    std::byte ReadScaledByte(size_t index, std::byte min_v, std::byte max_v) const override { return static_cast<std::byte>(static_cast<int>(min_v) + static_cast<int>(static_cast<double>(data[index])/alphabet_size*(static_cast<int>(max_v)-static_cast<int>(min_v)))); }
    bool ReadScaledBit(size_t index, bool min_v=0, bool max_v=1) const override { return (data[index] | min_v) & max_v; }
    emp::vector<int> MultiReadScaledInt(size_t index, size_t dist, int min_v, int max_v) const override {
      emp::vector<int> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledInt(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::vector<double> MultiReadScaledDouble(size_t index, size_t dist, double min_v, double max_v) const override {
      emp::vector<double> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledDouble(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::vector<std::byte> MultiReadScaledByte(size_t index, size_t dist, std::byte min_v, std::byte max_v) const override {
      emp::vector<std::byte> outvec(dist);
      for (size_t i= 0; i < dist; i ++) {
        outvec[i] = IsValid(index + i) ? ReadScaledByte(index + i, min_v, max_v) : min_v; 
      }
      return outvec;
    }
    emp::BitVector MultiReadScaledBit(size_t index, size_t dist, bool min_v=0, bool max_v=1) const override {
      if (min_v) {
        emp::BitVector outvec(dist, true);
        return outvec;
      }
      if (!max_v) {
        emp::BitVector outvec(dist, false);
        return outvec;
      }
      return data.Export(dist, index);
    }
    void WriteScaledInt(size_t index, int value, int min_v, int max_v) override { data[index] = static_cast<int>(static_cast<double>(value-min_v)*alphabet_size/(max_v-min_v)); }
    void WriteScaledDouble(size_t index, double value, double min_v, double max_v) override { data[index] = static_cast<int>((value-min_v)*alphabet_size/(max_v-min_v)); }
    void WriteScaledByte(size_t index, std::byte value, std::byte min_v, std::byte max_v) override { data[index] = static_cast<int>(static_cast<double>(static_cast<int>(value)-static_cast<int>(min_v))*alphabet_size/(static_cast<int>(max_v)-static_cast<int>(min_v))); }
    void WriteScaledBit(size_t index, bool value, bool min_v=0, bool max_v=1) override { data[index] = (value | min_v) & max_v;  }
    void MultiWriteScaledInt(size_t index, emp::vector<int> values, int min_v, int max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledInt(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledDouble(size_t index, emp::vector<double> values, double min_v, double max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledDouble(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledByte(size_t index, emp::vector<std::byte> values, std::byte min_v, std::byte max_v) override {
      size_t dist = values.size();
      for (size_t i= 0; i < dist; i ++) {
        if (IsValid(index + i)) WriteScaledByte(index + i, values[i], min_v, max_v); 
      }
    }
    void MultiWriteScaledBit(size_t index, emp::BitVector values, bool min_v=0, bool max_v=1) override {
      if (min_v) {
        values.SetAll();
      } else if (!max_v) {
        values.Clear();
      }
      data = data.Import( values, index );
    }
    
    
    
 



    size_t Mutate(emp::Random & random) override {
      
      mut_dist.Setup(mut_p, data.size()); // NEEDED HERE?????
      mut_sites.resize(data.size()); // NEEDED HERE?????

      const size_t num_muts = mut_dist.PickRandom(random);
      

      if (num_muts == 0) return 0;
      if (num_muts == 1) {
        const size_t pos = random.GetUInt(data.size());
        data.Toggle(pos);
        return 1;
      }

      // Only remaining option is num_muts > 1.
      mut_sites.Clear();
      for (size_t i = 0; i < num_muts; i++) {
        const size_t pos = random.GetUInt(data.size());
        if (mut_sites[pos]) { --i; continue; }  // Duplicate position; try again.
        mut_sites.Set(pos);
      }
      data ^= mut_sites;

      return num_muts;
    }

    emp::BitVector GetAllBits() const {return data;}


    size_t Mutate(emp::Random & random, double _mut) {
      mut_p = _mut;
      return Mutate(random);
    };

    void Clear() override { 
      data.Clear();
    }

    void SetAll() override { 
      data.SetAll();
    }

  };

  

  
  

}

#endif
