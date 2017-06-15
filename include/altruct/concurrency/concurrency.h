#pragma once

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

namespace altruct {
namespace concurrency {

#define LOCK(_mutex) \
    for (bool first = true; first;) \
        for (std::lock_guard<std::mutex> lock(_mutex); first; first = false)
            // { /* locked */  }

} // concurrency
} // altruct
