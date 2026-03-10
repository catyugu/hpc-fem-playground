#ifndef MPFEM_LOGGER_HPP
#define MPFEM_LOGGER_HPP

#include <chrono>
#include <mutex>
#include <string>
#include <sstream>
#include <cstdlib>

namespace mpfem {

/**
 * @brief Logging severity levels used by Logger.
 */
enum class LogLevel {
    Debug,
    Info,
    Warning,
    Error
};

/**
 * @brief Minimal thread-safe logger used by mpfem modules.
 */
class Logger {
public:
    /**
     * @brief Sets minimum severity threshold for output.
     */
    static void setLevel(LogLevel level);

    /**
     * @brief Logs a message when level is above configured threshold.
     */
    static void log(LogLevel level, const std::string &message);

    /**
     * @brief Returns elapsed milliseconds since program start.
     */
    static std::chrono::milliseconds elapsedMillis();

    /**
     * @brief Formats elapsed time as human-readable string.
     */
    static std::string formatElapsed();

private:
    static std::mutex mutex_;
    static LogLevel minimumLevel_;
    static std::chrono::steady_clock::time_point startTime_;
};

/**
 * @brief RAII timer that logs elapsed time on destruction.
 */
class ScopedTimer {
public:
    explicit ScopedTimer(const std::string &label, LogLevel level = LogLevel::Info);
    ~ScopedTimer();

    // Non-copyable, non-movable
    ScopedTimer(const ScopedTimer &) = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;
    ScopedTimer(ScopedTimer &&) = delete;
    ScopedTimer &operator=(ScopedTimer &&) = delete;

    void stop();
    double getElapsedSeconds() const;
    std::string getElapsedStr() const;

private:
    std::string label_;
    LogLevel level_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point end_;
    bool stopped_ = false;
};

/**
 * @brief Assertion macro that logs error and exits on failure.
 * 
 * Usage: Check(condition, "Error message");
 * If condition is false, logs error and calls std::exit(1).
 */
#define Check(cond, msg) \
    do { \
        if (!(cond)) { \
            ::mpfem::Logger::log(::mpfem::LogLevel::Error, msg); \
            std::exit(1); \
        } \
    } while (0)

} // namespace mpfem

#endif