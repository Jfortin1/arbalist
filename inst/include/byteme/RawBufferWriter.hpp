#ifndef BYTEME_RAW_BUFFER_WRITER_HPP
#define BYTEME_RAW_BUFFER_WRITER_HPP

#include "Writer.hpp"
#include <vector>

/**
 * @file RawBufferWriter.hpp
 *
 * @brief Write bytes to a raw buffer without any extra transformations.
 */

namespace byteme {

/**
 * @brief Write bytes to a raw buffer.
 *
 * This class will append bytes to an internal instance of a `std::vector` without any further transformations.
 * Not much else to say here.
 */
class RawBufferWriter : public Writer {
public:
    /**
     * @param reserve Initial size of the output buffer to reserve.
     */
    RawBufferWriter(size_t reserve = 0) {
        output.reserve(reserve);
    }

public:
    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        output.insert(output.end(), buffer, buffer + n);
    }

    void finish() {}

public:
    /**
     * @cond
     */
    // Exposed for back-compatibility only.
    std::vector<unsigned char> output;
    /**
     * @endcond
     */

    /**
     * @return Contents of the output buffer.
     * This should only be used after `finish()` is called.
     */
    std::vector<unsigned char>& get_output() {
        return output;
    }
};

}

#endif
