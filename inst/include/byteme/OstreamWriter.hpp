#ifndef BYTEME_OSTREAM_WRITER_HPP
#define BYTEME_OSTREAM_WRITER_HPP

#include <ostream>
#include <stdexcept>
#include "Writer.hpp"

/**
 * @file OstreamWriter.hpp
 *
 * @brief Write bytes to an arbitrary output stream.
 */

namespace byteme {

/**
 * @brief Read bytes from a `std::ostream`.
 *
 * @tparam Pointer_ A (possibly smart) pointer to an `std::ostream` object.
 *
 * This is just a wrapper around `std::ostream::write` for compatibility.
 */
template<class Pointer_ = std::ostream*>
class OstreamWriter : public Writer {
public:
    /**
     * @param output Pointer to an output stream.
     * This is assumed to live until `finish()` is called.
     */
    OstreamWriter(Pointer_ output) : my_output(std::move(output)) {}

    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        my_output->write(reinterpret_cast<const char*>(buffer), n);
        if (!(my_output->good())) {
            throw std::runtime_error("failed to write to arbitrary output stream");
        }
    }

    void finish() {
        my_output->flush();
        if (my_output->fail() || my_output->bad()) {
            throw std::runtime_error("failed to flush to arbitrary output stream");
        }
    }

private:
    Pointer_ my_output;
};

}

#endif
