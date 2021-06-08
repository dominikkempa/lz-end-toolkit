/**
 * @file    src/uint56.hpp
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

#ifndef __UINT56_HPP_INCLUDED
#define __UINT56_HPP_INCLUDED

#include <cstdint>
#include <limits>


class uint56 {
  private:
    std::uint32_t low;
    std::uint16_t mid;
    std::uint8_t high;

  public:
    uint56() {}
    uint56(std::uint32_t l, std::uint16_t m, std::uint8_t h)
      : low(l), mid(m), high(h) {}
    uint56(const uint56& a) : low(a.low), mid(a.mid), high(a.high) {}
    uint56(const std::int32_t& a) : low(a), mid(0), high(0) {}
    uint56(const std::uint32_t& a) : low(a), mid(0), high(0) {}
    uint56(const std::uint64_t& a) :
      low(a & 0xFFFFFFFF),
      mid((a >> 32) & 0xFFFF),
      high((a >> 48) & 0xFF) {}
    uint56(const std::int64_t& a) :
      low(a & 0xFFFFFFFFL),
      mid((a >> 32) & 0xFFFF),
      high((a >> 48) & 0xFF) {}

    inline operator uint64_t() const {
      return
        (((std::uint64_t)high) << 48) |
        (((std::uint64_t)mid) << 32) |
        (std::uint64_t)low; }
    inline bool operator == (const uint56& b) const {
      return (low == b.low) && (mid == b.mid) && (high == b.high); }
    inline bool operator != (const uint56& b) const {
      return (low != b.low) || (mid != b.mid) || (high != b.high); }
} __attribute__((packed));

namespace std {

template<>
struct is_unsigned<uint56> {
  public:
    static const bool value = true;
};

template<>
class numeric_limits<uint56> {
  public:
    static uint56 min() {
      return uint56(
          std::numeric_limits<std::uint32_t>::min(),
          std::numeric_limits<std::uint16_t>::min(),
          std::numeric_limits<std::uint8_t>::min());
    }

    static uint56 max() {
      return uint56(
          std::numeric_limits<std::uint32_t>::max(),
          std::numeric_limits<std::uint16_t>::max(),
          std::numeric_limits<std::uint8_t>::max());
    }
};

}  // namespace std

#endif  // __SRC_UINT56_HPP_INCLUDED
