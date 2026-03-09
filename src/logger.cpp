#include "logger.hpp"

#include <iomanip>
#include <iostream>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mpfem {

std::mutex Logger::mutex_;
LogLevel Logger::minimumLevel_ = LogLevel::Info;
std::chrono::steady_clock::time_point Logger::startTime_ = std::chrono::steady_clock::now();

void Logger::setLevel(LogLevel level)
{
    std::lock_guard<std::mutex> lock(mutex_);
    minimumLevel_ = level;
}

void Logger::log(LogLevel level, const std::string &message)
{
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (level < minimumLevel_) {
        return;
    }

    const char* ompiRank = std::getenv("OMPI_COMM_WORLD_RANK");
    if (ompiRank != nullptr && std::atoi(ompiRank) != 0) {
        return;
    }

    const char* pmiRank = std::getenv("PMI_RANK");
    if (pmiRank != nullptr && std::atoi(pmiRank) != 0) {
        return;
    }

#ifdef MFEM_USE_MPI
    int isInitialized = 0;
    MPI_Initialized(&isInitialized);
    if (isInitialized != 0) {
        int isFinalized = 0;
        MPI_Finalized(&isFinalized);
        if (isFinalized != 0) {
            return;
        }

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank != 0) {
            return;
        }
    }
#endif

    const char* levelStr = "";
    switch (level) {
        case LogLevel::Debug:   levelStr = "DEBUG"; break;
        case LogLevel::Info:    levelStr = "INFO"; break;
        case LogLevel::Warning: levelStr = "WARN"; break;
        case LogLevel::Error:   levelStr = "ERROR"; break;
    }

    std::cerr << "[" << levelStr << "] " << message << std::endl;
}

std::chrono::milliseconds Logger::elapsedMillis()
{
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime_);
}

std::string Logger::formatElapsed()
{
    auto ms = elapsedMillis().count();
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << (ms / 1000.0) << "s";
    return oss.str();
}

ScopedTimer::ScopedTimer(const std::string &label, LogLevel level)
    : label_(label)
    , level_(level)
    , start_(std::chrono::steady_clock::now())
{
}

ScopedTimer::~ScopedTimer()
{
    if (!stopped_) {
        stop();
    }
}

void ScopedTimer::stop()
{
    if (stopped_) return;
    stopped_ = true;
    end_ = std::chrono::steady_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_);
    double seconds = elapsed.count() / 1000000.0;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << label_ << " completed in " << seconds << "s";
    Logger::log(level_, oss.str());
}

double ScopedTimer::getElapsedSeconds() const
{
    auto endTime = stopped_ ? end_ : std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endTime - start_);
    return elapsed.count() / 1000000.0;
}

std::string ScopedTimer::getElapsedStr() const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << getElapsedSeconds() << "s";
    return oss.str();
}

} // namespace mpfem