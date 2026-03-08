#include "logger.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>

namespace mpfem {

namespace {

const char *const LOG_LEVEL_DEBUG = "DEBUG";
const char *const LOG_LEVEL_INFO = "INFO";
const char *const LOG_LEVEL_WARNING = "WARNING";
const char *const LOG_LEVEL_ERROR = "ERROR";

int logLevelToPriority(const LogLevel level)
{
    if (level == LogLevel::Debug) {
        return 0;
    }
    if (level == LogLevel::Info) {
        return 1;
    }
    if (level == LogLevel::Warning) {
        return 2;
    }
    return 3;
}

const char *toLevelText(const LogLevel level)
{
    if (level == LogLevel::Debug) {
        return LOG_LEVEL_DEBUG;
    }
    if (level == LogLevel::Info) {
        return LOG_LEVEL_INFO;
    }
    if (level == LogLevel::Warning) {
        return LOG_LEVEL_WARNING;
    }
    return LOG_LEVEL_ERROR;
}

} // namespace

std::mutex Logger::mutex_;
LogLevel Logger::minimumLevel_ = LogLevel::Info;
std::chrono::steady_clock::time_point Logger::startTime_ = std::chrono::steady_clock::now();

void Logger::setLevel(const LogLevel level)
{
    std::lock_guard<std::mutex> lock(mutex_);
    minimumLevel_ = level;
}

void Logger::log(const LogLevel level, const std::string &message)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (logLevelToPriority(level) < logLevelToPriority(minimumLevel_)) {
        return;
    }

    std::clog << "[" << toLevelText(level) << "] " << message << '\n';
}

void Logger::logWithTimestamp(const LogLevel level, const std::string &message)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (logLevelToPriority(level) < logLevelToPriority(minimumLevel_)) {
        return;
    }

    std::clog << "[" << toLevelText(level) << "] "
              << "[" << formatElapsed() << "] "
              << message << '\n';
}

std::chrono::milliseconds Logger::elapsedMillis()
{
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime_);
}

std::string Logger::formatElapsed()
{
    const auto ms = elapsedMillis().count();
    const auto secs = ms / 1000;
    const auto millis = ms % 1000;

    std::ostringstream oss;
    oss << std::setw(5) << secs << "." << std::setfill('0') << std::setw(3) << millis << "s";
    return oss.str();
}

ScopedTimer::ScopedTimer(const std::string &label, const LogLevel level)
    : label_(label), level_(level), start_(std::chrono::steady_clock::now())
{
}

ScopedTimer::~ScopedTimer()
{
    const auto end = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_);
    const double seconds = elapsed.count() / 1000.0;

    std::ostringstream oss;
    oss << label_ << " completed in " << std::fixed << std::setprecision(3) << seconds << "s";
    Logger::log(level_, oss.str());
}

} // namespace mpfem
