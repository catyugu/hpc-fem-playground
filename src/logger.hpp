#ifndef MPFEM_LOGGER_HPP
#define MPFEM_LOGGER_HPP

#include <chrono>
#include <mutex>
#include <string>

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
     * @param level New minimum level.
     */
    static void setLevel(LogLevel level);

    /**
     * @brief Logs a message when level is above configured threshold.
     * @param level Severity of current message.
     * @param message Text payload.
     */
    static void log(LogLevel level, const std::string &message);

    /**
     * @brief Logs a message with elapsed time since program start.
     * @param level Severity of current message.
     * @param message Text payload.
     */
    static void logWithTimestamp(LogLevel level, const std::string &message);

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

private:
    std::string label_;
    LogLevel level_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace mpfem

#endif
