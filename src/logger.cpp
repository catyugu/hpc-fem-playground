#include "logger.hpp"

#include <iostream>

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

} // namespace mpfem
