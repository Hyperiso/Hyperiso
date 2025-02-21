#ifndef HYPERISO_INTERFACE_H
#define HYPERISO_INTERFACE_H

#include "MemoryManager.h"
#include "WilsonInterface.h"
#include <thread>
#include <vector>
#include <functional>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <map>

class HyperisoInterface {
private:
    struct ThreadContext {
        std::unique_ptr<MemoryManager> memoryManager;
        std::unique_ptr<WilsonInterface> wilsonInterface;
    };

    static std::mutex instanceMutex;
    static std::map<std::thread::id, ThreadContext> threadInstances;

    HyperisoInterface() = default;

    static ThreadContext& getThreadContext(const std::string& lhaFile, const std::string& model, const std::vector<int>& models) {
        auto threadId = std::this_thread::get_id();

        std::lock_guard<std::mutex> lock(instanceMutex);

        if (threadInstances.find(threadId) == threadInstances.end()) {
            threadInstances[threadId] = ThreadContext{
                std::make_unique<MemoryManager>(lhaFile, models),
                std::make_unique<WilsonInterface>(model)
            };
        }

        return threadInstances[threadId];
    }

    class ThreadPool {
    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;

        std::mutex queueMutex;
        std::condition_variable condition;
        bool stop;

    public:
        explicit ThreadPool(size_t threads) : stop(false) {
            for (size_t i = 0; i < threads; ++i) {
                workers.emplace_back([this] {
                    while (true) {
                        std::function<void()> task;

                        {
                            std::unique_lock<std::mutex> lock(this->queueMutex);
                            this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });

                            if (this->stop && this->tasks.empty())
                                return;

                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }

                        task();
                    }
                });
            }
        }

        void enqueue(std::function<void()> task) {
            {
                std::lock_guard<std::mutex> lock(queueMutex);
                tasks.emplace(std::move(task));
            }
            condition.notify_one();
        }

        ~ThreadPool() {
            {
                std::lock_guard<std::mutex> lock(queueMutex);
                stop = true;
            }
            condition.notify_all();

            for (std::thread& worker : workers) {
                if (worker.joinable()) {
                    worker.join();
                }
            }
        }
    };

public:
    static HyperisoInterface& GetInstance() {
        static HyperisoInterface instance;
        return instance;
    }

    void Initialize(const std::string& lhaFile, const std::string& model, const std::vector<int>& models = {0}) {
        getThreadContext(lhaFile, model, models);
    }

    MemoryManager* GetMemoryManager() {
        return getThreadContext("", "", {}).memoryManager.get();
    }

    WilsonInterface* GetWilsonInterface() {
        return getThreadContext("", "", {}).wilsonInterface.get();
    }

    void RunInParallel(const std::vector<std::string>& lhaFiles, 
                    std::function<void(HyperisoInterface&, const std::string&)> task,
                    const std::string& model = "SM", 
                    const std::vector<int>& models = {0}, 
                    size_t maxThreads = std::thread::hardware_concurrency()) {
        ThreadPool threadPool(maxThreads);

        for (const auto& lhaFile : lhaFiles) {
            threadPool.enqueue([&, lhaFile]() {
                this->Initialize(lhaFile, model, models);
                task(*this, lhaFile);
            });
        }

    }

    void ClearThreadContext() {
        std::lock_guard<std::mutex> lock(instanceMutex);
        threadInstances.erase(std::this_thread::get_id());
    }

    HyperisoInterface(const HyperisoInterface&) = delete;
    HyperisoInterface& operator=(const HyperisoInterface&) = delete;
};

std::mutex HyperisoInterface::instanceMutex;
std::map<std::thread::id, HyperisoInterface::ThreadContext> HyperisoInterface::threadInstances;

#endif // HYPERISO_INTERFACE_H
