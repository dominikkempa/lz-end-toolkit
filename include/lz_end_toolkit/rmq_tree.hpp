/**
 * @file    rmq_tree.hpp
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

#ifndef __LZ_END_TOOLKIT_RMQ_TREE_HPP_INCLUDED
#define __LZ_END_TOOLKIT_RMQ_TREE_HPP_INCLUDED

#include <cstdint>
#include <limits>
#include <memory>
#include <algorithm>

#include "../utils/utils.hpp"


namespace lz_end_toolkit_private {

template<typename ValueType, int cache_line_size>
struct rmq_tree {
  public:
    typedef ValueType value_type;
    static_assert(sizeof(value_type) <= cache_line_size &&
        ((sizeof(value_type) & (sizeof(value_type) - 1))) == 0,
        "sizeof(value_type) in rmq_tree must be a power of two that "
        "does not exceed the cache line size");
    static_assert(
        cache_line_size > 0 &&
        (cache_line_size & (cache_line_size - 1)) == 0,
        "cache line size must be a positive integer that is a power of two");

  private:

    // The RMQ is implemented as a perfect balanced tree having the blk_size
    // degree. The tree in located in the m_data memory block. Each node of
    // the tree contains blk_size integer entries.
    //
    // The i-th leaf of the tree has blk_size entries that contain the
    // minimums of the blocks: m_tab[(i*blk_size + j)*blk_size ..
    // (i*blk_size + j + 1)*blk_size - 1] for j = 0,1,...,blk_size-1.
    // Last blocks might be truncated according to the size parameter.
    // If a block is empty, then the corresponding entry of the corresponding
    // leaf is equal to std::numeric_limits<value_type>::max().
    //
    // Each internal node of the tree also has blk_size entries. The j-th
    // entry contains the minimum of the blk_size entries of the corresponding
    // child node (which can be a leaf).
    //
    // Organization of the tree and navigation are implemented as follows.
    // The i-th leaf occupies the region m_data[m_offset + i*blk_size ..
    // m_offset + (i+1)*blk_size - 1]. So, the variable m_offset is used to
    // provide fast access to the leaves. If m_data[i] is an entry of a
    // non-root node, then m_data[i/blk_size-1] is an entry of the parent of
    // the node and it contains the minimum of all entries of the node
    // containing m_data[i]. Thus, the entries m_data[i & ~blk_mask] and
    // m_data[i | blk_mask] are, respectively, the first and the last entry
    // of the node containing the entry m_data[i]. (blk_mask is defined below)

    // blk_size is guaranteed to be a power of two.
    const std::uint64_t blk_size = cache_line_size / sizeof(value_type);
    const std::uint64_t blk_mask = blk_size - 1;

    const value_type *m_tab;	// the input array
    value_type *m_data;			  // the memory for the perfect balanced tree
    std::uint64_t m_offset;		// m_data[m_offset ..] contains the leaves

  public:
    rmq_tree(const value_type *tab, std::uint64_t size) {
      m_tab = tab;
      m_offset = 0;
      std::uint64_t leaves = blk_size;
      for (; leaves * blk_size <= size; leaves *= blk_size)
        m_offset += leaves;

      m_data =
        utils::aligned_allocate_array<value_type>(m_offset + leaves, 64);

      for (std::uint64_t i = 0; i < leaves; ++i) {
        value_type found = std::numeric_limits<value_type>::max();
        for (std::uint64_t j = i * blk_size;
            j < size && j < (i + 1) * blk_size; ++j)
          found = std::min(found, m_tab[j]);

        // (i % blk_size)-th entry of (i / blk_size)-th leaf
        m_data[m_offset + i] = found;
      }

      for (std::uint64_t i = m_offset + leaves; i > blk_size; i -= blk_size) {
        value_type found = std::numeric_limits<value_type>::max();
        for (std::uint64_t j = i - blk_size; j < i; ++j)
          found = std::min(found, m_data[j]);
        m_data[(i - 1) / blk_size - 1] = found;
      }
    }

    ~rmq_tree() {
      utils::aligned_deallocate(m_data);
    }

    // Return the boolean value telling whether the
    // minimum in the range is >= than given threshold.
    inline bool query(std::uint64_t beg,
        std::uint64_t end, std::uint64_t threshold) const {

      if (beg >= end) return false;

      std::uint64_t beg_blk_last = beg | blk_mask;
      std::uint64_t end_blk_first = end & ~blk_mask;

      if (beg_blk_last > end_blk_first)  // beg, end are inside the same block
        return !min_is_under_threshold(m_tab, beg, end, threshold);
      if (min_is_under_threshold(m_tab, beg, beg_blk_last + 1, threshold)
          || min_is_under_threshold(m_tab, end_blk_first, end, threshold))
        return false;
      if (beg_blk_last + 1 == end_blk_first)
        return true;

      // m_data[left] is entry with the minimum of
      //   m_tab[beg_blk_last-blk_size+1..beg_blk_last]
      // m_data[right] is entry with the minimum of
      //   m_tab[end_blk_first..end_blk_first+blk_size-1]
      std::uint64_t left = m_offset + beg / blk_size;
      std::uint64_t right = m_offset + end / blk_size;

      while ((left ^ right) & ~blk_mask) {

        // m_data[left] and m_data[right] ain't in the same node
        if (min_is_under_threshold(m_data,
              left + 1, (left | blk_mask) + 1, threshold))
          return false;
        if (min_is_under_threshold(m_data,
              right & ~blk_mask, right, threshold))
          return false;
        left = left / blk_size - 1;    // get parent
        right = right / blk_size - 1;  // get parent
      }
      return !min_is_under_threshold(m_data, left + 1, right, threshold);
    }

  private:
    inline static bool min_is_under_threshold(
        const value_type *array, std::uint64_t beg,
        std::uint64_t end, std::uint64_t threshold) {
      for (std::uint64_t i = beg; i < end; ++i)
        if ((std::uint64_t)array[i] < threshold)
          return true;
      return false;
    }
};

}  // namespace lz_end_toolkit_private

#endif  // __LZ_END_TOOLKIT_RMQ_TREE_HPP_INCLUDED
