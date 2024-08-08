#ifndef FRAGMENT_PARSER_HPP
#define FRAGMENT_PARSER_HPP

#include "byteme/SomeFileReader.hpp"
#include <vector>

template<class FragmentHandler>
class FragmentParser {
public:
    FragmentParser(FragmentHandler* h) : handler(h) {}

private:
    FragmentHandler* handler;

    std::string seq_name;
    int start_pos = 0;
    int end_pos = 0;
    std::string cell_name;
    int read_count = 0;

    unsigned char field = 0;
    bool empty = true; 
    size_t line_count = 1;

    void add_and_clear() {
        if (field != 4 || empty) {
            throw std::runtime_error("expected 5 non-empty fields (" + std::to_string(field + 1) + " detected) on line " + std::to_string(line_count));
        }
        handler->add(seq_name, start_pos, end_pos, cell_name, read_count, line_count);

        // Setting everything back.
        field = 0;
        empty = true;
        ++line_count;

        seq_name.clear();
        start_pos = 0;
        end_pos = 0;
        read_count = 0;
        cell_name.clear();
    }

    void parse_buffer(const char* ptr, size_t available) {
        for (size_t i = 0; i < available; ++i) {
            switch(ptr[i]) {
                case '\n':
                    add_and_clear();
                    break;

                case '\t':
                    if (empty) {
                        throw std::runtime_error("expected non-empty field " + std::to_string(field + 1) + " on line " + std::to_string(line_count));
                    }
                    ++field;
                    empty = true;
                    break;

                default:
                    switch (field) {
                        case 0:
                            seq_name += ptr[i];
                            break;
                        case 1:
                            if (!std::isdigit(ptr[i])) {
                                throw std::runtime_error("only digits should be present in the start position field on line " + std::to_string(line_count));
                            }
                            start_pos *= 10;
                            start_pos += ptr[i] - '0';
                            break;
                        case 2:
                            if (!std::isdigit(ptr[i])) {
                                throw std::runtime_error("only digits should be present in the end position field on line " + std::to_string(line_count));
                            }
                            end_pos *= 10;
                            end_pos += ptr[i] - '0';
                            break;
                        case 3:
                            cell_name += ptr[i];
                            break;
                    }
                    empty = false;
                    break;
            }
        }
    }

public:
    void run(const std::string& fragment_file) {
        byteme::SomeFileReader gzreader(fragment_file);
        bool remaining = false;

        /*** Reading through all the comments at the start. ***/
        bool in_comment = false;
        do {
            remaining = gzreader.load();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();
            auto ptr = reinterpret_cast<const char*>(buffer);

            bool breakout = false;
            for (size_t position = 0; position < available; ++position) {
                if (!in_comment) {
                    if (ptr[position] == '#') {
                        in_comment = true;
                    } else {
                        parse_buffer(ptr + position, available - position);
                        breakout = true;
                        break;
                    }
                } else {
                    if (ptr[position] == '\n') {
                        in_comment = false;
                        ++line_count;
                    }
                }
            }

            if (breakout) {
                break;
            }
        } while (remaining);

        /*** Now actually doing the parsing.  ***/
        do {
            remaining = gzreader.load();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();
            auto ptr = reinterpret_cast<const char*>(buffer);
            parse_buffer(ptr, available);
        } while (remaining);

        if (seq_name.size()) { // add the last record manually, when there's no trailing newline.
            add_and_clear();
        }
    }
};

#endif
