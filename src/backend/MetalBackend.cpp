#include "MetalBackend.h"

#include <algorithm>
#include <fstream>
#include <sstream>

#include "../inc/Metal/Metal/MTLTypes.hpp"

namespace Sibernetic {

namespace {
std::string readFile(const std::string &path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}
} // namespace

MetalBackend::MetalBackend(const char *libraryPath) {
  autoreleasePool_ = NS::AutoreleasePool::alloc()->init();

  device_ = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
  if (!device_) {
    if (autoreleasePool_) {
      autoreleasePool_->release();
      autoreleasePool_ = nullptr;
    }
    throw std::runtime_error("No Metal device found");
  }

  library_ = compileLibrary(libraryPath);

  queue_ = NS::TransferPtr(device_->newCommandQueue());
  if (!queue_) {
    throw std::runtime_error("Failed to create Metal command queue");
  }
}

MetalBackend::~MetalBackend() {
  pipelines_.clear();
  queue_.reset();
  library_.reset();
  device_.reset();
  if (autoreleasePool_) {
    autoreleasePool_->release();
  }
}

NS::SharedPtr<MTL::Library>
MetalBackend::compileLibrary(const std::string &sourcePath) {
  NS::Error *error = nullptr;
  std::string source = readFile(sourcePath);
  NS::String *nsSource =
      NS::String::string(source.c_str(), NS::UTF8StringEncoding);

  MTL::Library *lib = device_->newLibrary(nsSource, nullptr, &error);
  if (!lib) {
    std::string msg = "Failed to compile Metal library";
    if (error && error->localizedDescription()) {
      msg += ": ";
      msg += error->localizedDescription()->utf8String();
    }
    throw std::runtime_error(msg);
  }
  return NS::TransferPtr(lib);
}

MTL::ComputePipelineState *MetalBackend::getPipeline(const char *kernelName) {
  auto it = pipelines_.find(kernelName);
  if (it != pipelines_.end()) {
    return it->second.get();
  }

  auto pipeline = createPipeline(kernelName);
  auto *ptr = pipeline.get();
  pipelines_[kernelName] = std::move(pipeline);
  return ptr;
}

NS::SharedPtr<MTL::ComputePipelineState>
MetalBackend::createPipeline(const char *functionName) {
  MTL::Function *func = library_->newFunction(
      NS::String::string(functionName, NS::UTF8StringEncoding));
  if (!func) {
    throw std::runtime_error(std::string("Metal function not found: ") +
                             functionName);
  }
  auto funcPtr = NS::TransferPtr(func);

  NS::Error *error = nullptr;
  MTL::ComputePipelineState *pipeline =
      device_->newComputePipelineState(funcPtr.get(), &error);
  if (!pipeline) {
    std::string msg =
        "Failed to create pipeline for " + std::string(functionName);
    if (error && error->localizedDescription()) {
      msg += ": ";
      msg += error->localizedDescription()->utf8String();
    }
    throw std::runtime_error(msg);
  }
  return NS::TransferPtr(pipeline);
}

void MetalBackend::dispatch(
    const char *kernelName, uint32_t threadCount,
    std::function<void(MTL::ComputeCommandEncoder *)> setupArgs,
    bool waitForCompletion) {
  dispatch(getPipeline(kernelName), threadCount, std::move(setupArgs),
           waitForCompletion);
}

void MetalBackend::dispatch(
    MTL::ComputePipelineState *pipeline, uint32_t threadCount,
    std::function<void(MTL::ComputeCommandEncoder *)> setupArgs,
    bool waitForCompletion) {
  MTL::CommandBuffer *cmdBuf = queue_->commandBuffer();
  if (!cmdBuf) {
    throw std::runtime_error("Failed to create command buffer");
  }

  MTL::ComputeCommandEncoder *encoder = cmdBuf->computeCommandEncoder();
  if (!encoder) {
    throw std::runtime_error("Failed to create compute encoder");
  }

  encoder->setComputePipelineState(pipeline);
  setupArgs(encoder);

  NS::UInteger threads = static_cast<NS::UInteger>(threadCount);
  NS::UInteger tgWidth = std::max<NS::UInteger>(
      1, std::min<NS::UInteger>(threads,
                                pipeline->maxTotalThreadsPerThreadgroup()));

  encoder->dispatchThreads(MTL::Size::Make(threads, 1, 1),
                           MTL::Size::Make(tgWidth, 1, 1));
  encoder->endEncoding();

  cmdBuf->commit();
  if (waitForCompletion) {
    cmdBuf->waitUntilCompleted();
    if (cmdBuf->status() != MTL::CommandBufferStatusCompleted) {
      throw std::runtime_error(
          "Metal command buffer did not complete successfully");
    }
  }
}

void MetalBackend::finish() {
  // Metal command queue guarantees in-order execution, so submitting
  // an empty command buffer and waiting on it drains all prior work.
  MTL::CommandBuffer *cmdBuf = queue_->commandBuffer();
  if (!cmdBuf) {
    throw std::runtime_error("Failed to create command buffer for finish()");
  }
  cmdBuf->commit();
  cmdBuf->waitUntilCompleted();
  if (cmdBuf->status() != MTL::CommandBufferStatusCompleted) {
    std::string msg = "Metal GPU error during async execution";
    if (cmdBuf->error() && cmdBuf->error()->localizedDescription()) {
      msg += ": ";
      msg += cmdBuf->error()->localizedDescription()->utf8String();
    }
    throw std::runtime_error(msg);
  }
}

} // namespace Sibernetic
