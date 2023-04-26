#pragma once

#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <vector>

namespace altruct {
namespace concurrency {

/**
 * A synchronization helper.
 *
 *   LOCK(mu_mutex) {
 *     // do stuff
 *   }
 */
#define LOCK(_mutex) \
    for (bool first = true; first;) \
        for (std::lock_guard<std::mutex> lock(_mutex); first; first = false)
            // { /* locked */  }

/**
 * Parallelly executes all the jobs provided by the given job provider.
 * See `range_job_provider` for a job provider example.
 *
 * Jobs are exectued by the workers provided by the given worker provider.
 * A worker provider example:
 *       struct my_worker_provider {
 *           struct my_worker {
 *               // state (preallocated memory, precomputed tables, etc.)
 *
 *               JOB_RESULT execute_job(const JOB& job) {
 *                   // ...
 *               }
 *           };
 *           my_worker create_worker() { return my_worker(); }
 *       };
 *
 * Upon completion job result gets collected by the given result collector.
 * See `add_result_collector` for a result collector example.
 *
 * Each of the specified `num_threads` threads has exactly one worker.
 * It is assumed that obtaining jobs and collecting results is cheap
 * comparing to the actual job execution.
 *
 * Synchronization is taken care of in this method so none of the provided
 * objects need to worry about it. Thay do not have to be thread-safe.
 */
template<typename RESULT_COLLECTOR, typename JOB_PROVIDER, typename WORKER_PROVIDER>
void parallel_execute(RESULT_COLLECTOR& result_collector, JOB_PROVIDER& job_provider, WORKER_PROVIDER& worker_provider, int num_threads) {
    std::mutex job_provider_mutex;
    std::mutex result_collector_mutex;
    auto func = [&]() {
        auto worker = worker_provider.create_worker();
        while (true) {
            std::unique_lock<std::mutex> job_provider_lock(job_provider_mutex);
            if (!job_provider.has_next_job()) break;
            auto job = job_provider.next_job();
            job_provider_lock.unlock();
            auto job_result = worker.execute_job(job);
            std::unique_lock<std::mutex> result_collector_lock(result_collector_mutex);
            result_collector.collect_result(job_result, job);
        }
    };
    if (num_threads > 1) {
        std::vector<std::thread> t(num_threads);
        for (int i = 0; i < num_threads; ++i) {
            t[i] = std::thread(func);
        }
        for (int i = 0; i < num_threads; ++i) {
            t[i].join();
        }
    } else {
        func();
    }
}

/**
 * A result collector that simply accumulates values.
 */
template<typename T>
struct add_result_collector {
    T result = 0;
    template<typename JOB>
    void collect_result(const T& job_result, const JOB& job) { result += job_result; }
};

/**
 * A job provider that breaks range into smaller ones.
 *
 * E.g. `range_job_provider(12, 58, 10)` generates the following jobs:
 * [{12, 22}, {22, 32}, {32, 42}, {42, 52}, {52, 58}]
 */
template<typename I>
struct range_job_provider {
    I begin, end, len;
    range_job_provider(I begin, I end, I len) : begin(begin), end(end), len(len) {}
    bool has_next_job() { return begin < end; }
    std::pair<I, I> next_job() { begin += len; return{ begin - len, std::min(begin, end) }; }
};


} // concurrency
} // altruct
