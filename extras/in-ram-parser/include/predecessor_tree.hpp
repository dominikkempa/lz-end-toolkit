/**
 * @file    predecessor_tree.hpp
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

#ifndef __PREDECESSOR_TREE_HPP_INCLUDED
#define __PREDECESSOR_TREE_HPP_INCLUDED

#include <algorithm>
#include <cstdint>

#include "utils.hpp"


struct predecessor_tree {
  private:
    std::uint64_t m_size;      // size of the range
    std::uint64_t m_power;     // number of internal nodes at deepest level
    std::uint64_t m_node_cnt;  // number of internal nodes inc. root (id == 1)
    std::uint64_t m_word_cnt;  // number of of allocated words
    std::uint64_t *m_data;

  public:
    predecessor_tree(std::uint64_t size) {
      m_size = size;
      m_power = (1 << 6);
      m_node_cnt = m_power + 1;

      // Note: the number of internal
      // nodes is always at least 65.
      while ((m_power << 6) < m_size) {
        m_power <<= 6;
        m_node_cnt += m_power;
      }

      m_word_cnt =
        ((m_node_cnt - 1) >> 6) +        // internal nodes
        ((m_size + (1 << 6) - 1) >> 6);  // leaves

      // Aligned alloc guarantees address that
      // is aligned with cache line boundary.
      m_data = utils::aligned_allocate_array<std::uint64_t>(m_word_cnt, 64);
      std::fill(m_data, m_data + m_word_cnt, (std::uint64_t)0);
    }

    ~predecessor_tree() {
      utils::aligned_deallocate(m_data);
    }

  private:

    // Let b_0...b_63 be bits of the given word w where B_i
    // is (i+1)-th least significant bit. For a given i,
    // the function returns the largest j < i such that
    // B_j = 1, or 64 if there is no such j.
    inline static std::uint64_t bitpred(std::uint64_t w, std::uint64_t i) {
      if (i == 0) return 64;
      w &= (((std::uint64_t)1 << i) - 1);
      if (w == 0) return 64;
      else return 63 - __builtin_clzll(w);
    }

    // Analogous to bitpred, but returns successor.
    inline static std::uint64_t bitsucc(std::uint64_t w, std::uint64_t i) {
      if (i == 63) return 64;
      w &= (~(((std::uint64_t)1 << (i + 1)) - 1));
      if (w == 0) return 64;
      else return __builtin_ctzll(w);
    }

    inline void propagate(std::uint64_t cur_node_id) {
      while (cur_node_id > 1) {
        std::uint64_t leftmost_child_node_id = (cur_node_id << 6) - 62;
        std::uint64_t leftmost_child_word_id =
          ((leftmost_child_node_id - 2) >> 6);

        std::uint64_t cur_node_word_id = ((cur_node_id - 2) >> 6);
        std::uint64_t cur_node_offset = ((cur_node_id - 2) & 63);

        std::uint8_t curbit =
          ((m_data[cur_node_word_id] &
            ((std::uint64_t)1 << cur_node_offset)) > 0);
        std::uint8_t newbit =
          (m_data[leftmost_child_word_id] > 0);  // logical or of all bits

        if (curbit == newbit) break;
        else {
          m_data[cur_node_word_id] ^= ((std::uint64_t)1 << cur_node_offset);
          cur_node_id = ((cur_node_id + 62) >> 6);
        }
      }
    }

  public:

    // Bit corresponding to node/leaf with id x = 1, 2, ... (root is 1)
    // is stored in m_data[(x - 2) >> 6]. This is to guarantee that bits
    // of nodes having the same parent are inside the same word. Note
    // that this means we don't store the bit corresponding to the root.
    // This is ok, since this bit is never used anyway. Node also that
    // indexing of the argument i of this function starts from 0.
    inline void set(std::uint64_t i) {
      std::uint64_t node_id = m_node_cnt + i + 1;
      std::uint64_t word_id = ((node_id - 2) >> 6);
      std::uint64_t word_offset = ((node_id - 2) & 63);
      m_data[word_id] |= ((std::uint64_t)1 << word_offset);
      propagate((node_id + 62) >> 6);
    }

    inline void reset(std::uint64_t i) {
      std::uint64_t node_id = m_node_cnt + i + 1;
      std::uint64_t word_id = ((node_id - 2) >> 6);
      std::uint64_t word_offset = ((node_id - 2) & 63);
      m_data[word_id] &= (~((std::uint64_t)1 << word_offset));
      propagate((node_id + 62) >> 6);
    }

    void clear() {
      std::fill(m_data, m_data + m_word_cnt, (std::uint64_t)0);
    }

    inline std::uint64_t pred(std::uint64_t i) const {
      std::uint64_t cur_node_id = m_node_cnt + i + 1;
      while (cur_node_id > 1) {
        std::uint64_t word_id = ((cur_node_id - 2) >> 6);
        std::uint64_t word_offset = ((cur_node_id - 2) & 63);
        std::uint64_t bitpred_result = bitpred(m_data[word_id], word_offset);
        if (bitpred_result == 64) cur_node_id = ((cur_node_id + 62) >> 6);
        else {
          cur_node_id -= word_offset - bitpred_result;
          break;
        }
      }

      if (cur_node_id == 1)
        return 0;

      while (cur_node_id <= m_node_cnt) {
        std::uint64_t leftmost_child_node_id = (cur_node_id << 6) - 62;
        std::uint64_t leftmost_child_word_id =
          ((leftmost_child_node_id - 2) >> 6);
        std::uint64_t child_id =
          __builtin_clzll(m_data[leftmost_child_word_id]);
        cur_node_id = leftmost_child_node_id + (63 - child_id);
      }

      return cur_node_id - m_node_cnt;
    }

    inline std::uint64_t succ(std::uint64_t i) const {
      std::uint64_t cur_node_id = m_node_cnt + i + 1;
      while (cur_node_id > 1) {
        std::uint64_t word_id = ((cur_node_id - 2) >> 6);
        std::uint64_t word_offset = ((cur_node_id - 2) & 63);
        std::uint64_t bitsucc_result = bitsucc(m_data[word_id], word_offset);
        if (bitsucc_result == 64) cur_node_id = ((cur_node_id + 62) >> 6);
        else {
          cur_node_id += bitsucc_result - word_offset;
          break;
        }
      }

      if (cur_node_id == 1)
        return m_size;

      while (cur_node_id <= m_node_cnt) {
        std::uint64_t leftmost_child_node_id = (cur_node_id << 6) - 62;
        std::uint64_t leftmost_child_word_id =
          ((leftmost_child_node_id - 2) >> 6);
        std::uint64_t child_id =
          __builtin_ctzll(m_data[leftmost_child_word_id]);
        cur_node_id = leftmost_child_node_id + child_id;
      }

      return cur_node_id - m_node_cnt - 1;
    }
};

#endif  // __PREDECESSOR_TREE_HPP_INCLUDED
