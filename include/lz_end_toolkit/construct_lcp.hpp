/**
 * @file    construct_lcp.hpp
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

#ifndef __LZ_END_TOOLKIT_CONSTRUCT_LCP_HPP_INCLUDED
#define __LZ_END_TOOLKIT_CONSTRUCT_LCP_HPP_INCLUDED

#include <cstdint>

#include "../utils/utils.hpp"


namespace lz_end_toolkit_private {

template<typename char_type,
  typename sa_int_type,
  typename lcp_int_type>
void construct_lcp(
    const char_type *text,
    std::uint64_t text_length,
    const sa_int_type *sa,
    lcp_int_type *lcp_array) {

  lcp_int_type *phi = lcp_array;
  lcp_int_type *plcp = utils::allocate_array<lcp_int_type>(text_length);

  std::uint64_t phi_undefined_index = sa[0];
  for (std::uint64_t i = 1; i < text_length; ++i) {
    std::uint64_t addr = sa[i];
    phi[addr] = sa[i - 1];
  }

  std::uint64_t lcp = 0;
  for (std::uint64_t i = 0; i < text_length; ++i) {
    if (i == phi_undefined_index) {
      plcp[i] = 0;
      continue;
    }

    std::uint64_t j = phi[i];
    while (i + lcp < text_length && j + lcp < text_length
        && text[i + lcp] == text[j + lcp]) ++lcp;

    plcp[i] = (lcp_int_type)lcp;
    if (lcp > 0)
      --lcp;
  }

  for (std::uint64_t i = 0; i < text_length; ++i) {
    std::uint64_t addr = sa[i];
    lcp_array[i] = plcp[addr];
  }

  utils::deallocate(plcp);
}

}  // namespace lz_end_toolkit_private

#endif  // __LZ_END_TOOLKIT_CONSTRUCT_LCP_HPP_INCLUDED
