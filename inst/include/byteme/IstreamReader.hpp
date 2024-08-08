#ifndef BYTEME_ISTREAM_READER_HPP
#define BYTEME_ISTREAM_READER_HPP

#include <istream>
#include <vector>
#include <stdexcept>
#include "Reader.hpp"

/**
 * @file IstreamReader.hpp
 *
 * @brief Read bytes from an input stream.
 */

namespace byteme {

/**
 * @brief Read bytes from a `std::istream`.
 *
 * @tparam Pointer_ A (possibly smart) pointer to an `std::istream` object.
 *
 * This is just a wrapper around `std::istream::read`,
 * mostly to avoid having to remember the correct way to check for end of file.
 */
template<class Pointer_ = std::istream*>
class IstreamReader : public Reader {
public:
    /**
     * @param input Pointer to an input stream.
     * @param buffer_size Size of the buffer to use for reading.
     */
    IstreamReader(Pointer_ input, size_t buffer_size = 65536) : my_input(std::move(input)), my_buffer(buffer_size) {}

    bool load() {
        if (!my_okay) {
            return false;
        }

        my_input->read(reinterpret_cast<char*>(my_buffer.data()), my_buffer.size());
        my_read = my_input->gcount();

        if (my_read < my_buffer.size()) {
            if (my_input->eof()) {
                my_okay = false;
            } else {
                throw std::runtime_error("failed to finish reading the input stream");
            }
        }

        return true;
    }

    const unsigned char* buffer() const {
        return my_buffer.data();
    }

    size_t available() const {
        return my_read;
    }

private:
    Pointer_ my_input;
    std::vector<unsigned char> my_buffer;
    size_t my_read = 0;
    bool my_okay = true;
};

}

#endif
