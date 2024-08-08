#ifndef BYTEME_PERBYTE_HPP
#define BYTEME_PERBYTE_HPP

/**
 * @file PerByte.hpp
 * @brief Perform byte-by-byte extraction.
 */

#include <thread>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <algorithm>
#include <exception>
#include <type_traits>

#include "Reader.hpp"

namespace byteme {

/**
 * @cond
 */
template<class Pointer_>
void skip_zero_buffers(Pointer_& reader, size_t& available) {
    available = 0;
    while (reader->load()) {
        available = reader->available(); // continue collecting bytes if a zero-length buffer is returned without load() returning false.
        if (available) {
            break;
        }
    } 

    // If available == 0 on return, then reader->load() must be false,
    // and there are no more bytes left in the source.
}
/**
 * @endcond
 */

/**
 * @brief Perform byte-by-byte extraction from a `Reader` source.
 *
 * @tparam Type_ Character type to return, usually `char` for text or `unsigned char` for binary.
 * @tparam Pointer_ A (possibly smart) non-`const` pointer to a `Reader`.
 *
 * This wraps a `Reader` so that developers can avoid the boilerplate of managing blocks of bytes,
 * when all they want is to iterate over the bytes of the input.
 */
template<typename Type_ = char, class Pointer_ = Reader*>
class PerByte {
private:
    const Type_* my_ptr = nullptr;
    size_t my_available = 0;
    size_t my_current = 0;
    size_t my_overall = 0;

    Pointer_ my_reader;

    void refill() {
        my_overall += my_available;
        skip_zero_buffers(my_reader, my_available);
        my_ptr = reinterpret_cast<const Type_*>(my_reader->buffer());
        my_current = 0;
    }

public:
    /**
     * @param reader A (possibly smart) pointer to a `Reader` object.
     */
    PerByte(Pointer_ reader) : my_reader(std::move(reader)) {
        refill();
    }

    /**
     * @return Whether this instance still has bytes to be read.
     */
    bool valid() const {
        return my_current < my_available;
    }

    /**
     * Advance to the next byte, possibly reading a new chunk from the supplied `Reader`.
     * This should only be called if `valid()` is `true`.
     *
     * @return Whether this instance still has bytes to be read, i.e., the output of `valid()` after advancing to the next byte.
     */
    bool advance() {
        ++my_current;
        if (my_current < my_available) {
            return true;
        }

        refill();
        return my_available > 0; // Check that we haven't reached the end of the reader.
    }

    /**
     * @return The current byte.
     *
     * This should only be called if `valid()` is `true`.
     */
    Type_ get() const {
        return my_ptr[my_current];
    }

    /**
     * @return The position of the current byte since the start of the input.
     */
    size_t position() const {
        return my_overall + my_current;
    }

    /**
     * Extract up to `number` bytes from the `Reader` source and store them in the `output`.
     * This is equivalent to (but more efficient than) calling `get()` and then `advance()` up to `number` times,
     * only iterating while the return value of `advance()` is still true.
     * The number of successful iterations is returned in the output as the first pair element,
     * while the return value of the final `advance()` is returned as the second pair element.
     *
     * @param number Number of bytes to extract.
     * @param[out] output Pointer to an output buffer of length `number`.
     * This is filled with up to `number` bytes from the source.
     *
     * @return Pair containing (1) the number of bytes that were successfully read into `output`,
     * and (2) whether there are any more bytes available in the source for future `get()` or `extract()` calls.
     */
    std::pair<size_t, bool> extract(size_t number, Type_* output) {
        size_t original = number;
        bool okay = true;

        while (1) {
            auto start = my_ptr + my_current;
            auto leftover = my_available - my_current;

            if (leftover > number) {
                my_current += number;
                number = 0;
                std::copy(start, my_ptr + my_current, output);
                break;

            } else {
                number -= leftover;
                std::copy(start, my_ptr + my_available, output);
                refill();

                okay = (my_available > 0);
                if (number == 0 || !okay) {
                    break;
                }
                output += leftover;
            }
        }

        return std::make_pair(original - number, okay);
    }
};

/**
 * @brief Perform parallelized byte-by-byte extraction from a `Reader` source.
 *
 * @tparam Type_ Character type to return, usually `char` for text or `unsigned char` for binary.
 * @tparam Pointer_ A (possibly smart) non-`const` pointer to a `Reader`.
 *
 * This is much like `PerByte` except that the `Reader`'s loading operation is called in a separate thread,
 * thus allowing the caller to parse the bytes of the current chunk in parallel.
 */
template<typename Type_ = char, class Pointer_ = Reader*>
class PerByteParallel {
private:
    size_t my_current = 0;
    size_t my_available = 0;
    size_t my_overall = 0;

    Pointer_ my_reader;

    std::vector<Type_> my_buffer;
    size_t my_next_available = 0;
    bool my_finished = false;

private:
    std::thread my_thread;
    std::exception_ptr my_thread_err = nullptr;
    std::mutex my_mut;
    std::condition_variable my_cv;
    bool my_ready_input, my_ready_output;

    void thread_loop() {
        while (!my_finished) {
            std::unique_lock lck(my_mut);
            my_cv.wait(lck, [&]() { return my_ready_input; });
            my_ready_input = false;

            if (my_finished) { // an explicit kill signal from the destructor, see below.
                break;
            }

            try {
                skip_zero_buffers(my_reader, my_next_available);
                my_finished = my_next_available == 0; // see the definition of skip_zero_buffers().
            } catch (...) {
                my_thread_err = std::current_exception();
                my_finished = true;
            }

            my_ready_output = true;
            lck.unlock();
            my_cv.notify_one();
        }
    }

private:
    void refill() {
        std::unique_lock lck(my_mut);
        my_cv.wait(lck, [&]() { return my_ready_output; });
        my_ready_output = false;

        if (my_thread_err) {
            std::rethrow_exception(my_thread_err);
        }

        my_overall += my_available;

        auto ptr = reinterpret_cast<const Type_*>(my_reader->buffer());
        my_available = my_next_available;
        my_buffer.resize(my_available);
        my_current = 0;

        std::copy(ptr, ptr + my_available, my_buffer.begin());
        my_ready_input = true;

        lck.unlock();
        if (!my_finished) {
            my_cv.notify_one();
        }
    }

public:
    /**
     * @copydoc PerByte::PerByte()
     */
    PerByteParallel(Pointer_ reader) : my_reader(std::move(reader)) {
        my_ready_input = false;
        my_thread = std::thread([&]() { thread_loop(); });

        skip_zero_buffers(my_reader, my_next_available);
        my_ready_output = true; // run the first iteration of refill().
        refill();
    }

    /**
     * @cond
     */
    ~PerByteParallel() {
        if (!my_finished) {
            std::unique_lock lck(my_mut);
            my_finished = true;
            my_ready_input = true;
            lck.unlock(); // releasing the lock so that the notified thread doesn't immediately block.
            my_cv.notify_one();
        }
        my_thread.join();
    }
    /**
     * @endcond
     */

    /**
     * @copydoc PerByte::valid()
     */
    bool valid() const {
        return my_current < my_available;
    }

    /**
     * @copydoc PerByte::advance()
     */
    bool advance() {
        ++my_current;
        if (my_current < my_available) {
            return true;
        }

        if (my_finished) {
            return false;
        }
        refill();

        return my_available > 0; // confirm there's actually bytes to extract in the next round.
    }

    /**
     * @copydoc PerByte::get()
     */
    Type_ get() const {
        return my_buffer[my_current];
    }

    /**
     * @copydoc PerByte::position()
     */
    size_t position() const {
        return my_overall + my_current;
    }

    /**
     * Extract up to `number` bytes from the `Reader` source and store them in the `output`.
     * This is equivalent to (but more efficient than) calling `get()` and then `advance()` up to `number` times,
     * only iterating while the return value of `advance()` is still true.
     * The number of successful iterations is returned in the output as the first pair element,
     * while the return value of the final `advance()` is returned as the second pair element.
     *
     * @param number Number of bytes to extract.
     * @param[out] output Pointer to an output buffer of length `number`.
     * This is filled with up to `number` bytes from the source.
     *
     * @return Pair containing (1) the number of bytes that were successfully read into `output`,
     * and (2) whether there are any more bytes available in the source for future `get()` or `extract()` calls.
     */
    std::pair<size_t, bool> extract(size_t number, Type_* output) {
        size_t original = number;
        bool okay = true;

        while (1) {
            auto start = my_buffer.data() + my_current;
            auto leftover = my_available - my_current;

            if (leftover > number) {
                my_current += number;
                number = 0;
                std::copy(start, my_buffer.data() + my_current, output);
                break;

            } else {
                number -= leftover;
                std::copy(start, my_buffer.data() + my_available, output);

                if (my_finished) {
                    my_current += leftover;
                    okay = false;
                    break;
                }
                refill();

                okay = (my_available > 0);
                if (number == 0 || !okay) {
                    break;
                }
                output += leftover;
            }
        }

        return std::make_pair(original - number, okay);
    }
};

}

#endif
