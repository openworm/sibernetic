#pragma once

#include <cstdint>

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

// ============ Metal helpers ============

inline void bindBuffer(MTL::ComputeCommandEncoder *enc, MTL::Buffer *buf,
                       uint32_t index) {
  enc->setBuffer(buf, 0, index);
}

template <typename T>
inline void bindScalar(MTL::ComputeCommandEncoder *enc, const T &value,
                       uint32_t index) {
  enc->setBytes(&value, sizeof(T), index);
}


} // namespace Sibernetic
