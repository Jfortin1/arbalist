#ifndef BYTEME_SELF_CLOSING_FILE_HPP
#define BYTEME_SELF_CLOSING_FILE_HPP

#include <cstdio>
#include <string>
#include <stdexcept>

namespace byteme {

struct SelfClosingFILE {
public:
    SelfClosingFILE(const char * path, const char* mode) : handle(std::fopen(path, mode)) {
        if (!handle) {
            throw std::runtime_error("failed to open file at '" + std::string(path) + "'");
        }
        return;
    }

    ~SelfClosingFILE() {
        if (handle) {
            // Don't check return value as we can't throw here anyway.
            std::fclose(handle);
        }
        return;
    }

public:
    SelfClosingFILE(SelfClosingFILE&& x) {
        handle = x.handle;
        x.handle = NULL;
    }

    SelfClosingFILE& operator=(SelfClosingFILE&& x) {
        handle = x.handle;
        x.handle = NULL;
        return *this;
    }

    // Delete the remaining constructors.
    SelfClosingFILE(const SelfClosingFILE&) = delete;
    SelfClosingFILE& operator=(const SelfClosingFILE&) = delete;

public:
    FILE* handle;
};

}

#endif
