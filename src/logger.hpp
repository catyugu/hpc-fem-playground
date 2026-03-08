#ifndef MPFEM_LOGGER_HPP
#define MPFEM_LOGGER_HPP

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

private:
    static std::mutex mutex_;
    static LogLevel minimumLevel_;
};

} // namespace mpfem

#endif
