#pragma once

#include <functional>
#include <string>
#include <unordered_map>

#include "../inc/Metal/Foundation/NSSharedPtr.hpp"
#include "../inc/Metal/Metal/MTLCommandBuffer.hpp"
#include "../inc/Metal/Metal/MTLCommandQueue.hpp"
#include "../inc/Metal/Metal/MTLComputeCommandEncoder.hpp"
#include "../inc/Metal/Metal/MTLComputePipeline.hpp"
#include "../inc/Metal/Metal/MTLDevice.hpp"
#include "../inc/Metal/Metal/MTLLibrary.hpp"

namespace Sibernetic {

class MetalBackend {
public:
  static constexpr const char *kDefaultLibraryPath = "src/metal/sphFluid.metal";

  explicit MetalBackend(const char *libraryPath = kDefaultLibraryPath);
  ~MetalBackend();

  MetalBackend(const MetalBackend &) = delete;
  MetalBackend &operator=(const MetalBackend &) = delete;

  MTL::Device *device() const { return device_.get(); }
  MTL::CommandQueue *queue() const { return queue_.get(); }

  MTL::ComputePipelineState *getPipeline(const char *kernelName);

  void dispatch(const char *kernelName, uint32_t threadCount,
                std::function<void(MTL::ComputeCommandEncoder *)> setupArgs,
                bool waitForCompletion = true);
  void dispatch(MTL::ComputePipelineState *pipeline, uint32_t threadCount,
                std::function<void(MTL::ComputeCommandEncoder *)> setupArgs,
                bool waitForCompletion = true);

  /// Wait for all previously submitted command buffers to finish.
  void finish();

private:
  NS::SharedPtr<MTL::Library> compileLibrary(const std::string &sourcePath);
  NS::SharedPtr<MTL::ComputePipelineState>
  createPipeline(const char *functionName);

  NS::AutoreleasePool *autoreleasePool_ = nullptr;
  NS::SharedPtr<MTL::Device> device_;
  NS::SharedPtr<MTL::Library> library_;
  NS::SharedPtr<MTL::CommandQueue> queue_;
  std::unordered_map<std::string, NS::SharedPtr<MTL::ComputePipelineState>>
      pipelines_;
};

} // namespace Sibernetic
