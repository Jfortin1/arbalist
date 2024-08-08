#ifndef BYTEME_RAW_FILE_READER_HPP
#define BYTEME_RAW_FILE_READER_HPP

#include <vector>
#include <stdexcept>
#include <string>
#include <cstdio>
#include "Reader.hpp"
#include "SelfClosingFILE.hpp"

/**
 * @file RawFileReader.hpp
 *
 * @brief Read a file without any extra transformations.
 */

namespace byteme {

/**
 * @brief Read bytes from a file, usually text.
 *
 * This is basically a simple wrapper around `FILE` structures, with correct closing and error checking.
 * Mostly provided because I always forget how to interact with `ifstream` objects when I want a sequence of bytes.
 */
class RawFileReader : public Reader {
public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    RawFileReader(const char* path, size_t buffer_size = 65536) : my_file(path, "rb"), my_buffer(buffer_size) {}

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    RawFileReader(const std::string& path, size_t buffer_size = 65536) : RawFileReader(path.c_str(), buffer_size) {}

public:
    bool load() {
        if (!my_okay) {
            return false;
        }

        auto& handle = my_file.handle;
        my_read = std::fread(my_buffer.data(), sizeof(unsigned char), my_buffer.size(), handle);

        if (my_read < my_buffer.size()) {
            if (std::feof(handle)) {
                my_okay = false;
            } else {
                throw std::runtime_error("failed to read raw binary file (fread error " + std::to_string(std::ferror(handle)) + ")");
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
    SelfClosingFILE my_file;
    std::vector<unsigned char> my_buffer;
    size_t my_read = 0;
    bool my_okay = true;
};

}

#endif
