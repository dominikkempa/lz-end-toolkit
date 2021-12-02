/**
 * @file    main.cpp
 * @section LICENCE
 *
 * This file is part of LZ-End Toolkit v0.1.0
 * See: https://github.com/dominikkempa/lz-end-toolkit
 *
 * Published in:
 *   Dominik Kempa and Dmitry Kosolobov:
 *   LZ-End Parsing in Compressed Space.
 *   Data Compression Conference (DCC), IEEE, 2017.
 *
 * Copyright (C) 2016-2021
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <ctime>
#include <string>
#include <unistd.h>

#include "../../../include/compute_lzend.hpp"
#include "../../../include/space_efficient_vector.hpp"
#include "../../../include/types/uint40.hpp"
#include "../../../include/types/uint48.hpp"


template<typename char_type,
  typename text_offset_type>
void naive_parse(
    const char_type * const text,
    const std::uint64_t text_length,
    std::vector<lzend_phrase<char_type, text_offset_type, text_offset_type> > &parsing) {

  typedef lzend_phrase<char_type, text_offset_type, text_offset_type>
    phrase_type;

  std::vector<bool> phrase_end(text_length, false);
  std::uint64_t beg = 0;
  while (beg < text_length) {
    std::uint64_t source_link = 0;
    std::uint64_t max_phrase_len = 0;

    for (std::uint64_t phrase_len = 1;
        beg + phrase_len + 1 <= text_length &&
        phrase_len <= beg;
        ++phrase_len) {

      for (std::uint64_t i = 0; i <= beg - phrase_len; ++i) {
        if (phrase_end[i + phrase_len - 1]) {
          std::uint64_t lcp = 0;
          while (lcp < phrase_len && text[beg + lcp] == text[i + lcp])
            ++lcp;
          if (lcp == phrase_len) {
            source_link = 0;
            for (std::uint64_t tt = 0; tt < i + phrase_len - 1; ++tt)
              if (phrase_end[tt]) ++source_link;
            max_phrase_len = phrase_len;
          }
        }
      }
    }

    parsing.push_back(phrase_type());
    parsing.back().m_char = text[beg + max_phrase_len];
    parsing.back().m_len = max_phrase_len + 1;
    parsing.back().m_link = source_link;
    phrase_end[beg + max_phrase_len] = true;
    beg += max_phrase_len + 1;
  }
}

template<typename char_type,
  typename text_offset_type = std::uint32_t>
void test(
    char_type * const text,
    const std::uint64_t text_length) {

  // Declare types.
  typedef lzend_phrase<char_type, text_offset_type, text_offset_type> phrase_type;

  // Compute correct parsing.
  std::vector<phrase_type> parsing_correct;
  naive_parse(text, text_length, parsing_correct);

  // Compute parsing using tested algorithm.
  space_efficient_vector<phrase_type> parsing;
  compute_lzend<char_type, text_offset_type>(text, text_length, &parsing);
  std::uint64_t computed_parsing_size = parsing.size();

  // Compare sizes.
  if (parsing_correct.size() != computed_parsing_size) {
    fprintf(stderr, "\nError!\n");
    fprintf(stderr, "  text = ");
    for (std::uint64_t i = 0; i < text_length; ++i)
      fprintf(stderr, "%c", (char)text[i]);
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }
}

int main() {
  srand(time(0) + getpid());

  typedef std::uint8_t char_type;

#if 0

  // Test string over alphabet of size two.
  static const std::uint64_t max_text_length = 18;
  char_type *text = new char_type[max_text_length];
  for (std::uint64_t text_length = 1;
      text_length <= max_text_length; ++text_length) {
    fprintf(stderr, "text_length = %lu\n", text_length);
    std::uint64_t text_count = (1UL << text_length);
    for (std::uint64_t code = 0; code < text_count; ++code) {
      for (std::uint64_t j = 0; j < text_length; ++j)
        if (code & (1UL << j)) text[j] = 'a';
        else text[j] = 'b';

      test(text, text_length);
    }
  }
#else

  // Test string over alphabet of size three.
#ifdef NDEBUG
  static const std::uint64_t max_text_length = 13;
#else
  static const std::uint64_t max_text_length = 8;
#endif
  char_type *text = new char_type[max_text_length];
  for (std::uint64_t text_length = 1;
      text_length <= max_text_length; ++text_length) {
    fprintf(stderr, "text_length = %lu\n", text_length);
    std::uint64_t text_count = 1;
    for (std::uint64_t i = 0; i < text_length; ++i)
      text_count *= 3;
    for (std::uint64_t code = 0; code < text_count; ++code) {
      std::uint64_t x = code;
      for (std::uint64_t j = 0; j < text_length; ++j) {
        text[j] = 'a' + x % 3;
        x /= 3;
      }

      test(text, text_length);
    }
  }
#endif

  delete[] text;
  fprintf(stderr, "All test passed.\n");
}

