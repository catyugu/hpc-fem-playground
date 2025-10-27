#ifndef MFEM_PLUGIN_BASE_HPP
#define MFEM_PLUGIN_BASE_HPP

#include "mfem.hpp"
#include <string>
#include <memory>
#include <chrono>

namespace mfem_plugins
{

    // Performance metrics structure
    struct PerformanceMetrics
    {
        double setup_time_ms = 0.0;
        double solve_time_ms = 0.0;
        double total_time_ms = 0.0;
        int iterations = 0;
        double residual = 0.0;
        size_t memory_bytes = 0;

        void reset()
        {
            setup_time_ms = 0.0;
            solve_time_ms = 0.0;
            total_time_ms = 0.0;
            iterations = 0;
            residual = 0.0;
            memory_bytes = 0;
        }

        void print(const std::string &name) const
        {
            std::cout << "=== " << name << " Performance ===" << std::endl;
            std::cout << "Setup time:   " << setup_time_ms << " ms" << std::endl;
            std::cout << "Solve time:   " << solve_time_ms << " ms" << std::endl;
            std::cout << "Total time:   " << total_time_ms << " ms" << std::endl;
            std::cout << "Iterations:   " << iterations << std::endl;
            std::cout << "Residual:     " << residual << std::endl;
            std::cout << "Memory:       " << memory_bytes / (1024.0 * 1024.0) << " MB" << std::endl;
        }
    };

    // Timer utility
    class Timer
    {
    private:
        std::chrono::high_resolution_clock::time_point start_time;
        bool running;

    public:
        Timer() : running(false) {}

        void start()
        {
            start_time = std::chrono::high_resolution_clock::now();
            running = true;
        }

        double stop()
        {
            if (!running)
                return 0.0;
            auto end_time = std::chrono::high_resolution_clock::now();
            running = false;
            std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
            return elapsed.count();
        }

        double elapsed() const
        {
            if (!running)
                return 0.0;
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> elapsed = now - start_time;
            return elapsed.count();
        }
    };

    // Base plugin interface
    class PluginBase
    {
    protected:
        std::string name_;
        PerformanceMetrics metrics_;

    public:
        PluginBase(const std::string &name) : name_(name) {}
        virtual ~PluginBase() {}

        const std::string &name() const { return name_; }
        const PerformanceMetrics &metrics() const { return metrics_; }
        void reset_metrics() { metrics_.reset(); }

        virtual void print_info() const
        {
            std::cout << "Plugin: " << name_ << std::endl;
        }
    };

} // namespace mfem_plugins

#endif // MFEM_PLUGIN_BASE_HPP
