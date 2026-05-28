// This translation unit provides the metal-cpp private symbol definitions.
// These macros must be defined before including metal-cpp headers in exactly
// one .cpp file per linked binary. See inc/Metal/README.md for details.
#ifdef SIBERNETIC_USE_METAL

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#include "../inc/Metal/Foundation/Foundation.hpp" // IWYU pragma: keep
#include "../inc/Metal/Metal/Metal.hpp"           // IWYU pragma: keep

#endif // SIBERNETIC_USE_METAL
