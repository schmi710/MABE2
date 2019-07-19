/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019
 *
 *  @file  ConfigType.h
 *  @brief Maintains information about a single type.
 *  @note Status: ALPHA
 */

#ifndef MABE_CONFIG_TYPE_H
#define MABE_CONFIG_TYPE_H

#include "base/assert.h"
#include "base/vector.h"

namespace mabe {

  enum class BaseType {
    VOID,
    VALUE,
    STRING,
    STRUCT
  };

  class ConfigType {
  private:
    BaseType base_type; ///< What is the base type?
    size_t array_depth; ///< How nested is the array?

    bool IsBaseType() const { return array_depth == 0; }
    bool IsArray() const { return array_depth > 0; }

    bool IsVoid() const { return IsBaseType() && base_type == BaseType::VOID; };
    bool IsValue() const { return IsBaseType() && base_type == BaseType::VOID; };
    bool IsString() const { return IsBaseType() && base_type == BaseType::VOID; };
    bool IsStruct() const { return IsBaseType() && base_type == BaseType::VOID; };
  };

}

#endif