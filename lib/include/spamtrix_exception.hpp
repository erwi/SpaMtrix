// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
#define SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
#include <stdexcept>
#include <string>

#define ERROR_LOCATION std::string(__FILE__) + ":" + std::to_string(__LINE__) + " in " + std::string(__func__)

/**
 * @brief Exception type used by SpaMtrix.
 */
class SpaMtrixException : public std::runtime_error {
public:
    /**
     * @brief Construct an exception with the supplied message.
     *
     * @param msg Error message.
     */
    explicit SpaMtrixException(const std::string &msg) : std::runtime_error(msg) {}
};

#endif //SPAMTRIX_SPAMTRIX_EXCEPTION_HPP
