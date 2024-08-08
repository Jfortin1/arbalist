#ifndef BYTEME_GZIP_FILE_READER_HPP
#define BYTEME_GZIP_FILE_READER_HPP

#include "zlib.h"
#include <stdexcept>
#include <vector>
#include <string>
#include "SelfClosingGzFile.hpp"
#include "Reader.hpp"

/**
 * @file GzipFileReader.hpp
 *
 * @brief Read a Gzip-compressed file.
 */

namespace byteme {

/**
 * @brief Read uncompressed bytes from a Gzip-compressed file.
 *
 * This is basically a wrapper around Zlib's `gzFile` with correct closing and error checking.
 */
class GzipFileReader : public Reader {
public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    GzipFileReader(const char* path, size_t buffer_size = 65536) : my_gzfile(path, "rb"), my_buffer(buffer_size) {}

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    GzipFileReader(const std::string& path, size_t buffer_size = 65536) : GzipFileReader(path.c_str(), buffer_size) {}

public:
    bool load() {
        my_read = gzread(my_gzfile.handle, my_buffer.data(), my_buffer.size());
        if (my_read) {
            return true;
        }

        if (!gzeof(my_gzfile.handle)) { 
            int dummy;
            throw std::runtime_error(gzerror(my_gzfile.handle, &dummy));
        }

        return false;
    }

    const unsigned char* buffer() const {
        return my_buffer.data();
    }

    size_t available() const {
        return my_read;
    }

private:
    SelfClosingGzFile my_gzfile;
    std::vector<unsigned char> my_buffer;
    size_t my_read = 0;
};

}

#endif
