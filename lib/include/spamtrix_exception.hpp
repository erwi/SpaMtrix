#ifndef SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
#define SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
#include <stdexcept>

#define ERROR_LOCATION std::string(__FILE__) + ":" + std::to_string(__LINE__) + " in " + std::string(__func__)

class SpaMtrixException : public std::runtime_error {
public:
    explicit SpaMtrixException(const std::string &msg) : std::runtime_error(msg) {}
};

#endif //SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
