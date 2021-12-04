/**
 * @file    compute_lzend.hpp
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

#ifndef __COMPUTE_LZEND_HPP_INCLUDED
#define __COMPUTE_LZEND_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "space_efficient_vector.hpp"
#include "predecessor_tree.hpp"
#include "rmq_tree.hpp"
#include "construct_lcp.hpp"
#include "construct_sa.hpp"
#include "utils.hpp"


template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
struct lzend_phrase {
  char_type m_char;          // last symbol of a phrase
  text_offset_type m_link;   // source phrase id
  phrase_offset_type m_len;  // length of phrase, including last symbol
} __attribute__((packed));

template<typename char_type,
  typename text_offset_type,
  typename lcp_int_type>
bool marked_lcp(
    const std::uint64_t text_pos,
    const std::uint64_t text_length,
    const std::uint64_t absorb_length,
    const std::uint64_t * const sa_of_rev_text,
    const text_offset_type * const isa_of_rev_text,
    const predecessor_tree &marked,
    const space_efficient_vector<text_offset_type> &recent_phrase_ends,
    std::uint64_t &prev_phrase_id,
    const rmq_tree<lcp_int_type, 64> &lcp_rmq) {

  const std::uint64_t pos = isa_of_rev_text[text_pos];
  std::uint64_t left = marked.pred(pos);

  if (left > 0) {
    if (lcp_rmq.query(left, pos + 1, absorb_length)) {
      std::uint64_t low = 0;
      std::uint64_t high = recent_phrase_ends.size();
      const std::uint64_t key =
        text_length - 1 - (std::uint64_t)sa_of_rev_text[left - 1];

      while (high - low > 1) {
        const std::uint64_t mid = (low + high) / 2;
        if ((std::uint64_t)recent_phrase_ends[mid] > key)
          high = mid;
        else low = mid;
      }

      prev_phrase_id = low;
      return true;
    }
  }

  std::uint64_t right = marked.succ(pos);
  if (right < text_length) {
    if (lcp_rmq.query(pos + 1, right + 1, absorb_length)) {
      std::uint64_t low = 0;
      std::uint64_t high = recent_phrase_ends.size();
      const std::uint64_t key =
        text_length - 1 - (std::uint64_t)sa_of_rev_text[right];

      while (high - low > 1) {
        const std::uint64_t mid = (low + high) / 2;
        if ((std::uint64_t)recent_phrase_ends[mid] > key)
          high = mid;
        else low = mid;
      }

      prev_phrase_id = low;
      return true;
    }
  }

  return false;
}

template<typename char_type,
  typename text_offset_type>
void compute_lzend(
    char_type * const text,
    const std::uint64_t text_length,
    const std::uint64_t * const sa_of_rev_text,
    space_efficient_vector<
      lzend_phrase<char_type, text_offset_type, text_offset_type> > *parsing,
    bool verbose = false) {

  // Allocate suffix array, inverse suffix array,
  // LCP array and predecessor tree.
  typedef std::uint64_t lcp_int_type;
  text_offset_type * const isa_of_rev_text =
    utils::allocate_array<text_offset_type>(text_length);
  lcp_int_type * const lcp_array =
    utils::allocate_array<lcp_int_type>(text_length);
  predecessor_tree * const marked = new predecessor_tree(text_length);

  // Process blocks left-to-right.
  typedef lzend_phrase<char_type, text_offset_type, text_offset_type>
    phrase_type;

  // Compute LCP array.
  {
    long double compute_lcp_start = (long double)0;
    if (verbose == true) {
      fprintf(stderr, "Compute LCP: ");
      compute_lcp_start = utils::wclock();
    }
    std::reverse(text, text + text_length);
    construct_lcp<char_type, std::uint64_t, lcp_int_type>(text,
        text_length, sa_of_rev_text, lcp_array);
    std::reverse(text, text + text_length);
    if (verbose == true) {
      long double elapsed = utils::wclock() - compute_lcp_start;
      fprintf(stderr, "%.2Lfs\n", elapsed);
    }
  }

  // Compute ISA array.
  {
    long double compute_isa_start = (long double)0;
    if (verbose == true) {
      fprintf(stderr, "Compute ISA: ");
      compute_isa_start = utils::wclock();
    }
    for (std::uint64_t i = 0; i < text_length; ++i) {
      const std::uint64_t addr =
        text_length - 1 - (std::uint64_t)sa_of_rev_text[i];
      isa_of_rev_text[addr] = i;
    }
    if (verbose == true) {
      long double elapsed = utils::wclock() - compute_isa_start;
      fprintf(stderr, "%.2Lfs\n", elapsed);
    }
  }

  // Initialize RMQ for LCP array.
  typedef rmq_tree<lcp_int_type, 64> rmq_type;
  space_efficient_vector<text_offset_type> recent_phrase_ends;
  rmq_type *lcp_rmq = NULL;
  {
    long double rmq_init_start = (long double)0;
    if (verbose == true) {
      fprintf(stderr, "Initialize RMQ: ");
      rmq_init_start = utils::wclock();
    }
    lcp_rmq = new rmq_type(lcp_array, text_length);
    if (verbose == true) {
      long double elapsed = utils::wclock() - rmq_init_start;
      fprintf(stderr, "%.2Lfs\n", elapsed);
    }
  }

  // Parse the text.
  long double parsing_start = (long double)0;
  if (verbose == true) {
    parsing_start = utils::wclock();
    fprintf(stderr, "Parse: ");
  }
  for (std::uint64_t parsed_prefix_length = 0;
      parsed_prefix_length < text_length; ++parsed_prefix_length) {

    // Extend the parsing to include text[i].
    std::uint64_t cur_parsing_size = parsing->size();
    std::uint64_t prev_phrase_id = 0;

    if (cur_parsing_size >= 1) {

      // Try absorbing the last two phrases.
      bool found = false;
      if (cur_parsing_size >= 2) {
        const std::uint64_t last_phrase_len =
          (*parsing)[cur_parsing_size - 1].m_len;
        const std::uint64_t two_last_phrases_len = last_phrase_len +
          (std::uint64_t)(*parsing)[cur_parsing_size - 2].m_len;

        const std::uint64_t text_pos_to_ignore =
          parsed_prefix_length - last_phrase_len - 1;
        const std::uint64_t isa_pos_to_ignore =
          isa_of_rev_text[text_pos_to_ignore];
        marked->reset(isa_pos_to_ignore);
        if (marked_lcp<char_type, text_offset_type, lcp_int_type>(
              parsed_prefix_length - 1, text_length,
              two_last_phrases_len, sa_of_rev_text,
              isa_of_rev_text, *marked, recent_phrase_ends, prev_phrase_id,
              *lcp_rmq)) {
          parsing->pop_back();
          recent_phrase_ends.pop_back();
          parsing->back().m_len = two_last_phrases_len + 1;
          found = true;
        } else {
          marked->set(isa_pos_to_ignore);
        }
      }

      // Try absorbing the last phrase.
      if (!found) {
        const std::uint64_t last_phrase_len =
          (*parsing)[cur_parsing_size - 1].m_len;
        if (marked_lcp<char_type, text_offset_type, lcp_int_type>(
              parsed_prefix_length - 1, text_length,
              last_phrase_len,
              sa_of_rev_text, isa_of_rev_text,
              *marked, recent_phrase_ends,
              prev_phrase_id, *lcp_rmq)) {
          parsing->back().m_len = last_phrase_len + 1;
          found = true;
        }
      }

      // Do not absorb any phrases.
      if (!found) {
        parsing->push_back(phrase_type());
        parsing->back().m_len = 1;
        recent_phrase_ends.push_back(text_offset_type());
      }

    } else {
      parsing->push_back(phrase_type());
      parsing->back().m_len = 1;
      recent_phrase_ends.push_back(text_offset_type());
    }


    if ((std::uint64_t)(parsing->back().m_len) == 1
        && parsing->size() > 1) {
      const std::uint64_t text_pos_to_set = parsed_prefix_length - 1;
      const std::uint64_t isa_pos_to_set = isa_of_rev_text[text_pos_to_set];
      marked->set(isa_pos_to_set);
    }

    parsing->back().m_char = text[parsed_prefix_length];
    parsing->back().m_link = prev_phrase_id;
    recent_phrase_ends.back() = parsed_prefix_length;
  }
  if (verbose == true) {
    long double parsing_time = utils::wclock() - parsing_start;
    fprintf(stderr, "%.2Lfs\n", parsing_time);
  }

  // Clean up.
  delete lcp_rmq;
  delete marked;
  utils::deallocate(lcp_array);
  utils::deallocate(isa_of_rev_text);
}

template<typename char_type,
  typename text_offset_type>
void compute_lzend(
    char_type * const text,
    const std::uint64_t text_length,
    space_efficient_vector<lzend_phrase<char_type, text_offset_type, text_offset_type> > *parsing,
    bool verbose = false) {

  // Allocate the suffix array of reversed text.
  std::uint64_t * const sa_of_rev_text =
    utils::allocate_array<std::uint64_t>(text_length);

  // Compute SA.
  {
    long double compute_sa_start = (long double)0;
    if (verbose == true) {
      fprintf(stderr, "Compute SA: ");
      compute_sa_start = utils::wclock();
    }
    std::reverse(text, text + text_length);
    construct_sa<char_type, std::uint64_t>(
        text, text_length, sa_of_rev_text);
    std::reverse(text, text + text_length);
    if (verbose == true) {
      long double elapsed = utils::wclock() - compute_sa_start;
      fprintf(stderr, "%.2Lfs\n", elapsed);
    }
  }

  // Compute the parsing.
  compute_lzend<char_type, text_offset_type>(text, text_length,
      sa_of_rev_text, parsing, verbose);

  // Clean up.
  utils::deallocate(sa_of_rev_text);
}

template<typename char_type,
  typename text_offset_type>
void compute_lzend_for_file(
    std::string text_filename,
    std::string output_filename,
    bool verbose = true) {

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  output_filename = utils::absolute_path(output_filename);

  // Initialize I/O stats and start the timer.
  long double parsing_start = utils::wclock();
  utils::initialize_stats();

  // Compute text length.
  std::uint64_t text_length = utils::file_size(text_filename);

  // Print initial message.
  if (verbose == true) {
    fprintf(stderr, "Compute LZ-End parsing\n");
    fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
    fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
    fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
    fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length,
        (1.L * text_length * sizeof(char_type)) / (1UL << 20));
    fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
    fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
        sizeof(text_offset_type));
    fprintf(stderr, "\n\n");
  }

  // Allocate and read text from file.
  char_type *text = utils::allocate_array<char_type>(text_length);
  utils::read_from_file(text, text_length, text_filename);

  // Compute the parsing.
  typedef lzend_phrase<char_type, text_offset_type, text_offset_type>
    phrase_type;
  typedef space_efficient_vector<phrase_type> vector_type;
  vector_type *parsing = new vector_type();
  compute_lzend<char_type, text_offset_type>(text,
      text_length, parsing, verbose);

  // Write parsing to file.
  {

    // Allocate buffer.
    static const std::uint64_t bufsize = (64 << 10);
    typedef lzend_phrase<char_type, text_offset_type, text_offset_type>
      output_phrase_type;
    output_phrase_type *buf =
      utils::allocate_array<output_phrase_type>(bufsize);

    // Open output file.
    std::FILE *f = utils::file_open(output_filename, "w");

    // Write the header.
    {
      std::uint64_t char_type_bits = 8 * sizeof(char_type);
      std::uint64_t text_offset_type_bits = 8 * sizeof(text_offset_type);
      std::uint64_t header = 0;
      header |= (char_type_bits - 1);
      header |= ((text_offset_type_bits - 1) << 8);
      utils::write_to_file(&header, 1, f);
    }

    // Write data to disk.
    std::uint64_t written = 0;
    while (written < parsing->size()) {
      std::uint64_t filled =
        std::min(bufsize, parsing->size() - written);
      for (std::uint64_t i = 0; i < filled; ++i) {
        buf[i].m_char = (*parsing)[written + i].m_char;
        buf[i].m_link = (*parsing)[written + i].m_link;
        buf[i].m_len = (std::uint64_t)((*parsing)[written + i].m_len);
      }
      written += filled;
      utils::write_to_file(buf, filled, f);
    }

    // Clean up.
    std::fclose(f);
    utils::deallocate(buf);
  }

  // Store number of phrases.
  const std::uint64_t n_phrases = parsing->size();

  // Compute length of the longest phrase.
  std::uint64_t max_phrase_length = 0;
  for (std::uint64_t i = 0; i < n_phrases; ++i) {
    const std::uint64_t len = (*parsing)[i].m_len;
    max_phrase_length = std::max(max_phrase_length, len);
  }

  // Clean up.
  delete parsing;
  utils::deallocate(text);

  // Print summary.
  if (verbose == true) {
    long double elapsed = utils::wclock() - parsing_start;
    long double rel_time = (elapsed * 1000000.L) / text_length;
    long double avg_phrase_length = (1.L * text_length) / n_phrases;
    fprintf(stderr, "\n\nComputation finished. Summary:\n");
    fprintf(stderr, "  Absolute time = %.2Lfs\n", elapsed);
    fprintf(stderr, "  Relative time = %.2Lfus/symbol\n", rel_time);
    fprintf(stderr, "  Number of phrases = %lu\n", n_phrases);
    fprintf(stderr, "  Maximal phrase length = %lu\n", max_phrase_length);
    fprintf(stderr, "  Average phrase length = %.2Lf\n", avg_phrase_length);
    fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
        utils::get_current_ram_allocation(),
        (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
  }
}

#endif  // __COMPUTE_LZEND_HPP_INCLUDED
