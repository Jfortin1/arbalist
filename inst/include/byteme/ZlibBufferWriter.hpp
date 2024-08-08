#ifndef BYTEME_ZLIB_BUFFER_WRITER_HPP
#define BYTEME_ZLIB_BUFFER_WRITER_HPP

#include "zlib.h"
#include <stdexcept>
#include <vector>
#include "Writer.hpp"

/**
 * @file ZlibBufferWriter.hpp
 *
 * @brief Write bytes to a Zlib-compressed buffer.
 */

namespace byteme {

/**
 * @brief Compress and write bytes to a Zlib-compressed buffer.
 *
 * This is basically a wrapper around Zlib's deflate method, with correct closing and error checking.
 */
class ZlibBufferWriter : public Writer {
private:
    /**
     * @cond
     */
    struct ZStream {
        ZStream(int mode, int level) {
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.avail_in = 0;
            strm.next_in = Z_NULL;

            // Check out https://zlib.net/manual.html for the constants.
            int windowBits;
            if (mode == 0) { // deflate
                windowBits = -15;
            } else if (mode == 1) { // Zlib
                windowBits = 15;
            } else if (mode == 2) { // Gzip
                windowBits = 16 + 15;
            } else {
                throw std::runtime_error("unknown Zlib compression mode supplied");
            } 

            int ret = deflateInit2(&strm, level, Z_DEFLATED, windowBits, 8, Z_DEFAULT_STRATEGY);
            if (ret != Z_OK) {
                throw std::runtime_error("failed to initialize Zlib buffer compression");
            }
        }

        ~ZStream() {
            deflateEnd(&strm);
            return;
        }

    public:
        // Delete the remaining constructors.
        ZStream(const ZStream&) = delete;
        ZStream(ZStream&&) = delete;
        ZStream& operator=(const ZStream&) = delete;
        ZStream& operator=(ZStream&&) = delete;

    public:
        z_stream strm;
    };
    /**
     * @endcond
     */

public:
    /**
     * @param mode Compression mod of the stream - DEFLATE (0), Zlib (1) or Gzip (2).
     * @param compression_level Compression level from 1 to 9.
     * Larger values improve compression at the cost of speed.
     * @param buffer_size Size of the compression buffer in bytes.
     * Larger values improve speed at the cost of memory.
     */
    ZlibBufferWriter(int mode = 2, int compression_level = 6, size_t buffer_size = 65536) : my_zstr(mode, compression_level), my_holding(buffer_size) {}

public:
    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        my_zstr.strm.next_in = const_cast<unsigned char*>(buffer); // for C compatibility.
        my_zstr.strm.avail_in = n;
        dump(Z_NO_FLUSH);
    }

    void finish() {
        my_zstr.strm.next_in = nullptr;
        my_zstr.strm.avail_in = 0;
        dump(Z_FINISH);
    }

private:
    ZStream my_zstr;
    std::vector<unsigned char> my_holding;

    void dump(int flag) {
        do {
            my_zstr.strm.avail_out = my_holding.size();
            my_zstr.strm.next_out = my_holding.data();
            deflate(&(my_zstr.strm), flag); // no need to check, see https://zlib.net/zlib_how.html.
            size_t compressed = my_holding.size() - my_zstr.strm.avail_out;
            output.insert(output.end(), my_holding.begin(), my_holding.begin() + compressed);
        } while (my_zstr.strm.avail_out == 0);
    }

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
