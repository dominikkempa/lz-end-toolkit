/**
 * @file    src/lz_end_toolkit_src/decode.hpp
 * @section LICENCE
 *
 * This file is part of LZ-End Toolkit v0.1.0
 * Published in:
 *   Dominik Kempa and Dmitry Kosolobov:
 *   LZ-End Parsing in Compressed Space.
 *   Data Compression Conference (DCC), IEEE, 2017.
 *
 * Copyright (C) 2016-2017
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *   Dmitry Kosolobov <dkosolobov (at) mail.ru>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __SRC_LZ_END_TOOLKIT_SRC_DECODE_HPP_INCLUDED
#define __SRC_LZ_END_TOOLKIT_SRC_DECODE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <stack>
#include <deque>
#include <string>
#include <algorithm>

#include "io/async_stream_writer.hpp"
#include "../uint24.hpp"
#include "../uint40.hpp"
#include "../uint48.hpp"
#include "../uint56.hpp"
#include "utils.hpp"


namespace lz_end_toolkit_private {

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
struct phrase {
  char_type m_char;
  text_offset_type m_link;
  phrase_offset_type m_len;
} __attribute__((packed));

template<typename text_offset_type>
struct stack_item {
  text_offset_type m_len;
  text_offset_type m_id;

  stack_item() {}
  stack_item(text_offset_type len, text_offset_type id)
    : m_len(len),
      m_id(id) {}
} __attribute__((packed));

template<typename char_type,
  typename text_offset_type>
void decode(
    std::string parsing_filename,
    std::string output_filename) {
  typedef phrase<char_type, text_offset_type,
    text_offset_type> phrase_type;

  // Start the timer.
  long double start = utils::wclock();

  // Turn paths absolute.
  parsing_filename = utils::absolute_path(parsing_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print basic info.
  fprintf(stderr, "Running LZ-End Toolkit v0.1.0\n");
  fprintf(stderr, "Mode = decoding\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Parsing filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "Output (text) filename = %s\n", output_filename.c_str());
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));

  // Read the parsing from disk.
  static const std::uint64_t header_size_bytes = 8;
  std::uint64_t parsing_file_size = utils::file_size(parsing_filename);
  std::uint64_t parsing_size =
    (parsing_file_size - header_size_bytes) / sizeof(phrase_type);
  phrase_type *parsing = new phrase_type[parsing_size];
  utils::read_at_offset(parsing,
      header_size_bytes, parsing_size, parsing_filename);

  // Compute text length.
  std::uint64_t text_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < parsing_size; ++phrase_id)
    text_length += (std::uint64_t)parsing[phrase_id].m_len;

  // Check if text length is nonzero.
  if (text_length == 0) {
    fprintf(stderr, "\nError: text length is 0\n");
    std::exit(EXIT_FAILURE);
  }

  // Check if a text this long could have been
  // encoded using the given integer type.
  if ((std::uint64_t)std::numeric_limits<text_offset_type>::max() <
      text_length - 1) {
    fprintf(stderr, "\nError: sizeof(text_offset_type) = %lu "
        "is too small.\n", sizeof(text_offset_type));
    std::exit(EXIT_FAILURE);
  }

  // Initialize the output (text) streamer.
  typedef async_stream_writer<char_type> output_writer_type;
  output_writer_type *output_writer =
    new output_writer_type(output_filename);

  // Print parsing size.
  fprintf(stderr, "Number of phrases = %lu\n", parsing_size);
  fprintf(stderr, "\n\n");

  // Decode the parsing and write to disk.
  std::uint64_t decoded_prefix_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < parsing_size; ++phrase_id) {
    if ((phrase_id % 10000) == 0)
      fprintf(stderr, "\rDecode: %.2Lf%%",
          (100.L * decoded_prefix_length) / text_length);

    // Check if phrase length is not zero.
    if ((std::uint64_t)parsing[phrase_id].m_len == 0) {
      fprintf(stderr, "\nError: phrase length has length zero!\n");
      std::exit(EXIT_FAILURE);
    }

    // Check if the phrase link is correct.
    if ((std::uint64_t)parsing[phrase_id].m_len > 1 &&
        (std::uint64_t)parsing[phrase_id].m_link >= phrase_id) {
      fprintf(stderr, "\nError: incorrect link of the phrase!\n");
      fprintf(stderr, "\n  Prev phrase link for phrase[%lu] is %lu\n",
          phrase_id, (std::uint64_t)parsing[phrase_id].m_link);
      std::exit(EXIT_FAILURE);
    }

    typedef stack_item<text_offset_type> stack_item_type;
    typedef std::deque<stack_item_type> deque_type;
    typedef std::stack<stack_item_type, deque_type> stack_type;

    // Decode next phrase.
    stack_type s;
    s.push(stack_item_type(
          (std::uint64_t)parsing[phrase_id].m_len,
          phrase_id));

    while (!s.empty()) {
      stack_item_type x = s.top();
      s.pop();

      if (x.m_len > 1) {
        s.push(stack_item_type((text_offset_type)1, x.m_id));
        std::uint64_t rest = x.m_len - 1;
        std::uint64_t id = x.m_id;
        if (parsing[x.m_id].m_len > 1) {
          s.push(stack_item_type(
                std::min(rest, (std::uint64_t)parsing[id].m_len - 1),
                parsing[id].m_link));

          rest -= std::min(rest, (std::uint64_t)parsing[id].m_len - 1);
        }
        --id;

        while (rest > 0) {
          s.push(stack_item_type(
                std::min(rest, (std::uint64_t)parsing[id].m_len),
                id));
          rest -= std::min(rest, (std::uint64_t)parsing[id].m_len);
          --id;
        }
      } else {
        char_type next_decoded_char = parsing[x.m_id].m_char;
        output_writer->write(next_decoded_char);
        ++decoded_prefix_length;
      }
    }
  }

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rDecode: 100.00%%\n\n\n");
  fprintf(stderr, "Decoding finished. Summary:\n");
  fprintf(stderr, "  Absolute time = %.2Lfs\n", total_time);
  fprintf(stderr, "  Decoded text length = %lu\n",
      decoded_prefix_length);

  // Clean up.
  delete output_writer;
  delete[] parsing;
}

template<typename char_type>
void decode(
    std::string parsing_filename,
    std::string output_filename,
    std::uint64_t text_offset_type_bits) {
  if (text_offset_type_bits == 32)
    decode<char_type, std::uint32_t>(parsing_filename, output_filename);
  else if (text_offset_type_bits == 40)
    decode<char_type, uint40>(parsing_filename, output_filename);
  else if (text_offset_type_bits == 48)
    decode<char_type, uint48>(parsing_filename, output_filename);
  else if (text_offset_type_bits == 56)
    decode<char_type, uint56>(parsing_filename, output_filename);
  else if (text_offset_type_bits == 64)
    decode<char_type, std::uint64_t>(parsing_filename, output_filename);
  else {
    fprintf(stderr, "\nError: decode: unsupported "
        "text_offset_type_bits (%lu)!\n", text_offset_type_bits);
    std::exit(EXIT_FAILURE);
  }
}

void decode(
    std::string parsing_filename,
    std::string output_filename,
    std::uint64_t char_type_bits,
    std::uint64_t text_offset_type_bits) {
  if (char_type_bits == 8)
    decode<std::uint8_t>(parsing_filename,
        output_filename, text_offset_type_bits);
  else {
    fprintf(stderr, "\nError: decode: unsupported "
        "char_type_bits (%lu)!\n", char_type_bits);
    std::exit(EXIT_FAILURE);
  }
}

void decode(
    std::string parsing_filename,
    std::string output_filename) {

  // Decode header.
  std::uint64_t char_type_bits = 0;
  std::uint64_t text_offset_type_bits = 0;
  {
    std::FILE *f = utils::file_open(parsing_filename, "r");
    std::uint64_t header = 0;
    utils::read_from_file(&header, 1, f);
    char_type_bits = (header & ((std::uint64_t)0xFF)) + 1;
    text_offset_type_bits = ((header >> 8) & ((std::uint64_t)0xFF)) + 1;
    std::fclose(f);
  }

  // Run the verification.
  decode(parsing_filename, output_filename,
      char_type_bits, text_offset_type_bits);
}

}  // namespace lz_end_toolkit_private

void decode(
    std::string parsing_filename,
    std::string output_filename) {
  lz_end_toolkit_private::decode(parsing_filename, output_filename);
}

#endif  // __SRC_LZ_END_TOOLKIT_SRC_DECODE_HPP_INCLUDED
