/**
 * @file    parse.hpp
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

#ifndef __LZ_END_TOOLKIT_PARSE_HPP_INCLUDED
#define __LZ_END_TOOLKIT_PARSE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <limits>
#include <algorithm>
#include <deque>
#include <stack>
#include <unistd.h>

#include "hash_table.hpp"
#include "space_efficient_vector.hpp"
#include "predecessor_tree.hpp"
#include "rmq_tree.hpp"
#include "construct_lcp.hpp"
#include "construct_sa.hpp"
#include "../utils/utils.hpp"
#include "../types/uint24.hpp"
#include "../types/uint40.hpp"
#include "../types/uint48.hpp"


namespace lz_end_toolkit_private {

template<typename S, typename T>
struct packed_pair {
  typedef packed_pair<S, T> pair_type;

  packed_pair() {}
  packed_pair(S &f, T &s) {
    first = f;
    second = s;
  }

  packed_pair(S f, T s) {
    first = f;
    second = s;
  }

  inline bool operator == (const pair_type &p) const {
    return first == p.first && second == p.second;
  }

  S first;
  T second;
} __attribute__((packed));

template<typename first_type, typename second_type>
std::uint64_t get_hash(const packed_pair<first_type, second_type> &x) {
  return (std::uint64_t)x.first * (std::uint64_t)29996224275833 +
    (std::uint64_t)x.second * (std::uint64_t)14638944639703;
}

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
struct phrase {
  char_type m_char;          // last symbol of a phrase
  text_offset_type m_link;   // source phrase id
  phrase_offset_type m_len;  // length of phrase, including last symbol
} __attribute__((packed));

template<typename char_type,
  typename text_offset_type,
  typename node_id_type>
struct trie_node {
  typedef trie_node<char_type,
    text_offset_type,
    node_id_type> node_type;

  char_type m_char;               // first symbol on the edge from the parent
  text_offset_type m_depth;       // string-length of path from the root
  text_offset_type m_phr;         // phrase ID corresponding to any leaf below
  node_id_type m_parent;          // parent of the node (NULL for root)
  node_id_type m_leftmost_child;  // leftmost child (can be NULL)
  node_id_type m_right_sibling;   // next child of the parent (can be NULL)

  trie_node() {
    m_parent = std::numeric_limits<node_id_type>::max();
    m_leftmost_child = std::numeric_limits<node_id_type>::max();
    m_right_sibling = std::numeric_limits<node_id_type>::max();
  }

  std::uint64_t get_child(
      char_type c,
      const space_efficient_vector<node_type> &trie_nodes) const {
    std::uint64_t child = m_leftmost_child;
    while (child != std::numeric_limits<node_id_type>::max()) {
      if (trie_nodes[child].m_char == c) return child;
      else child = trie_nodes[child].m_right_sibling;
    }

    return std::numeric_limits<std::uint64_t>::max();
  }

  inline void add_child(
      std::uint64_t newchild,
      space_efficient_vector<node_type> &trie_nodes) {
    if (m_leftmost_child == std::numeric_limits<node_id_type>::max())
      m_leftmost_child = newchild;
    else {
      std::uint64_t child = m_leftmost_child;
      while (trie_nodes[child].m_right_sibling !=
          std::numeric_limits<node_id_type>::max())
        child = trie_nodes[child].m_right_sibling;
      trie_nodes[child].m_right_sibling = newchild;
    }
  }

  inline void remove_child(
      std::uint64_t todelete,
      space_efficient_vector<node_type> &trie_nodes) {
    if ((std::uint64_t)m_leftmost_child == todelete)
      m_leftmost_child = trie_nodes[m_leftmost_child].m_right_sibling;
    else {
      std::uint64_t child = m_leftmost_child;
      while ((std::uint64_t)trie_nodes[child].m_right_sibling != todelete)
        child = trie_nodes[child].m_right_sibling;
      trie_nodes[child].m_right_sibling =
        trie_nodes[trie_nodes[child].m_right_sibling].m_right_sibling;
    }
    trie_nodes[todelete].m_right_sibling =
      std::numeric_limits<node_id_type>::max();
  }

  inline bool has_child() const {
    return m_leftmost_child != std::numeric_limits<node_id_type>::max();
  }

  inline bool has_one_child(
      const space_efficient_vector<node_type> &trie_nodes) const {

    return m_leftmost_child != std::numeric_limits<node_id_type>::max() &&
      trie_nodes[(std::uint64_t)m_leftmost_child].m_right_sibling ==
      std::numeric_limits<node_id_type>::max();
  }
} __attribute__((packed));

std::uint64_t rst(
    std::uint64_t x,
    std::uint64_t i) {
  return (x & (~((1UL << i) - 1)));
}

std::uint64_t compute_pv(
    std::uint64_t d,
    std::uint64_t pd) {
  std::uint64_t i = 0;
  while (rst(d, i + 1) > pd) ++i;
  return rst(d, i);
}

// Return (a * b) mod p, where p = (2^k) - 1.
// Requires a, b <= 2^k. Tested for k = 1, .., 63.
std::uint64_t mul_mod_meresenne(std::uint64_t a,
    std::uint64_t b, std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  unsigned __int128 ab =
    (unsigned __int128)a *
    (unsigned __int128)b;
  std::uint64_t lo = (std::uint64_t)ab;
  std::uint64_t hi = (ab >> 64);
  lo = (lo & p) + ((lo >> k) + (hi << (64 - k)));
  lo = (lo & p) + (lo >> k);
  return lo == p ? 0 : lo;
}

// Return a mod p, where p = (2^k) - 1.
// Works for any a in [0..2^64).
// Tested for k = 1, .., 63.
std::uint64_t mod_mersenne(std::uint64_t a, std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  if (k < 32) {

    // We need to check if a <= 2^(2k).
    std::uint64_t threshold = ((std::uint64_t)1 << (k << 1));
    if (a <= threshold) {
      a = (a & p) + (a >> k);
      a = (a & p) + (a >> k);
      return a == p ? 0 : a;
    } else return a % p;
  } else {

    // We are guaranteed that a < 2^(2k)
    // because a < 2^64 <= 2^(2k).
    a = (a & p) + (a >> k);
    a = (a & p) + (a >> k);
    return a == p ? 0 : a;
  }
}

// Return random number x in [0..p), where p = (2^k) - 1.
inline std::uint64_t rand_mod_mersenne(std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  return utils::random_int64((std::uint64_t)0, p - 1);
}

template<typename char_type>
std::uint64_t hash_of_reversed_substr(
    std::uint64_t beg,
    std::uint64_t end,
    std::uint64_t text_length,
    std::uint64_t hash_window_offset,
    std::uint64_t mersenne_prime_exponent,
    const std::uint64_t *hash_power,
    const std::uint64_t *hash_prefix_of_rev_text) {

  std::uint64_t mersenne_prime = 
    ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
  std::uint64_t a =
    hash_prefix_of_rev_text[text_length - beg - hash_window_offset];
  std::uint64_t b = hash_power[end - beg];
  std::uint64_t c =
    hash_prefix_of_rev_text[text_length - end - hash_window_offset];
  std::uint64_t bc_mod_p =
    mul_mod_meresenne(b, c, mersenne_prime_exponent);
  if (a >= bc_mod_p)
    return a - bc_mod_p;
  else
    return mersenne_prime - (bc_mod_p - a);
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type>
std::uint64_t approx_find(
    const char_type *text_window,
    std::uint64_t window_offset,
    std::uint64_t pat_end,
    std::uint64_t pat_length,
    std::uint64_t trie_root,
    std::uint64_t text_length,
    std::uint64_t hash_window_offset,
    std::uint64_t mersenne_prime_exponent,
    const std::uint64_t *hash_power,
    const std::uint64_t *hash_prefix_of_rev_text,
    const space_efficient_vector<trie_node<char_type,
        text_offset_type, node_id_type> > &trie_nodes,
    const hash_table<packed_pair<phrase_offset_type,
        std::uint64_t>, node_id_type, node_id_type> &nav) {

  typedef packed_pair<phrase_offset_type, std::uint64_t> nav_key_type;
  typedef node_id_type nav_value_type;

  std::uint64_t p = 0;
  std::uint64_t v = trie_root;

  std::uint64_t pat_length_logceil = 0;
  std::uint64_t pow2 = 1;
  while (pow2 < pat_length) {
    pow2 <<= 1;
    ++pat_length_logceil;
  }

  for (std::uint64_t iplus = pat_length_logceil + 1; iplus > 0; --iplus) {
    std::uint64_t i = iplus - 1;

    if ((std::uint64_t)trie_nodes[v].m_depth >= p + (1UL << i))
      p += (1UL << i);
    else if (p + (1UL << i) <= pat_length) {
      std::uint64_t pp = p + (1UL << i);
      std::uint64_t hh = hash_of_reversed_substr<char_type>(pat_end - pp,
          pat_end, text_length, hash_window_offset, mersenne_prime_exponent,
          hash_power, hash_prefix_of_rev_text);

      const nav_value_type *ret =
        nav.find(nav_key_type((phrase_offset_type)pp, hh));

      if (ret != NULL) {
        p += (1UL << i);
        v = *ret;
      }
    }
  }

  if ((std::uint64_t)trie_nodes[v].m_depth < pat_length) {
    std::uint64_t vnext = trie_nodes[v].get_child(text_window[pat_end - 1 -
        (std::uint64_t)trie_nodes[v].m_depth - window_offset], trie_nodes);
    if (vnext != std::numeric_limits<std::uint64_t>::max())
      v = vnext;
  }

  return trie_nodes[v].m_phr;
}

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
std::uint64_t longest_common_suffix(
    const char_type *pat_window,
    std::uint64_t window_offset,
    std::uint64_t pat_length,
    std::uint64_t max_lcs,
    std::uint64_t phrase_id,
    const space_efficient_vector<phrase<char_type,
        text_offset_type, phrase_offset_type> > &parsing,
    char_type &mismatch_char) {

  typedef std::pair<text_offset_type, text_offset_type> stack_item_type;
  typedef std::deque<stack_item_type> deque_type;
  typedef std::stack<stack_item_type, deque_type> stack_type;

  stack_type s;
  std::uint64_t lcs = 0;
  s.push(std::make_pair(
        (text_offset_type)phrase_id,
        (text_offset_type)pat_length));

  while (lcs < pat_length && lcs < max_lcs && !s.empty()) {
    stack_item_type x = s.top();
    std::uint64_t id = x.first;
    std::uint64_t len = x.second;
    s.pop();

    if (parsing[id].m_char ==
        pat_window[pat_length - 1 - lcs - window_offset]) {

      ++lcs;
      std::uint64_t rest = len - 1;
      if (rest > 0) {
        bool topush = false;
        stack_item_type item_to_push =
          std::make_pair((text_offset_type)0, (text_offset_type)0);
        if (parsing[id].m_len > 1) {
          std::uint64_t prev_len =
            std::min(rest, (std::uint64_t)parsing[id].m_len - 1);
          std::uint64_t prev_id = parsing[id].m_link;
          item_to_push = std::make_pair(
              (text_offset_type)prev_id,
              (text_offset_type)prev_len);
          topush = true;
          rest -= prev_len;
        }

        if (rest > 0 && id > 0) {
          --id;
          s.push(std::make_pair(
                (text_offset_type)id,
                (text_offset_type)rest));
        }

        if (topush)
          s.push(item_to_push);
      }
    } else {
      mismatch_char = parsing[id].m_char;
      break;
    }
  }

  return lcs;
}

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
std::uint64_t hash_of_reversed_substr(
    std::uint64_t phrase_id,
    std::uint64_t length,
    std::uint64_t mersenne_prime_exponent,
    std::uint64_t hash_variable,
    const space_efficient_vector<phrase<char_type,
        text_offset_type, phrase_offset_type> > &parsing) {

  typedef std::pair<text_offset_type, text_offset_type> stack_item_type;
  typedef std::deque<stack_item_type> deque_type;
  typedef std::stack<stack_item_type, deque_type> stack_type;

  stack_type s;
  std::uint64_t ret = 0;
  s.push(std::make_pair(
        (text_offset_type)phrase_id,
        (text_offset_type)length));

  while (!s.empty()) {
    stack_item_type x = s.top();
    std::uint64_t id = x.first;
    std::uint64_t len = x.second;
    s.pop();

    {
      std::uint64_t ab_mod_p = mul_mod_meresenne(
          ret, hash_variable, mersenne_prime_exponent);
      std::uint64_t c_mod_p = mod_mersenne(
          (std::uint64_t)parsing[id].m_char,
          mersenne_prime_exponent);
      ret = mod_mersenne(ab_mod_p + c_mod_p,
          mersenne_prime_exponent);
    }

    std::uint64_t rest = len - 1;
    if (rest > 0) {
      bool topush = false;
      stack_item_type item_to_push =
        std::make_pair((text_offset_type)0, (text_offset_type)0);

      if (parsing[id].m_len > 1) {
        std::uint64_t prev_len =
          std::min(rest, (std::uint64_t)parsing[id].m_len - 1);
        std::uint64_t prev_id = parsing[id].m_link;
        item_to_push = std::make_pair(
            (text_offset_type)prev_id,
            (text_offset_type)prev_len);
        topush = true;
        rest -= prev_len;
      }

      if (rest > 0 && id > 0) {
        --id;
        s.push(std::make_pair(
              (text_offset_type)id,
              (text_offset_type)rest));
      }

      if (topush)
        s.push(item_to_push);
    }
  }

  return ret;
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type>
std::uint64_t nearest_common_ancestor(
    std::uint64_t x,
    std::uint64_t y,
    const space_efficient_vector<trie_node<char_type,
        text_offset_type, node_id_type> > &trie_nodes) {

  while (x != y) {
    if ((std::uint64_t)trie_nodes[x].m_depth >=
        (std::uint64_t)trie_nodes[y].m_depth)
      x = (std::uint64_t)trie_nodes[x].m_parent;
    else y = (std::uint64_t)trie_nodes[y].m_parent;
  }

  return x;
}

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type,
  typename window_offset_type>
bool marked_lcp(
    std::uint64_t text_pos,
    std::uint64_t text_window_offset,
    std::uint64_t window_size,
    std::uint64_t absorb_length,
    std::uint64_t isa_pos_to_ignore,
    const window_offset_type *sa_of_rev_text_window,
    const window_offset_type *isa_of_rev_text_window,
    const predecessor_tree &marked,
    const space_efficient_vector<window_offset_type> &recent_phrase_ends,
    std::uint64_t nonrecent_phrase_count,
    std::uint64_t &prev_phrase_id,
    const rmq_tree<window_offset_type, 64> &lcp_rmq) {
  std::uint64_t pos = isa_of_rev_text_window[text_pos - text_window_offset];

  std::uint64_t left = marked.pred(pos);
  if (left > 0 && left == isa_pos_to_ignore + 1)
    left = marked.pred(isa_pos_to_ignore);
  if (left > 0) {
    if (lcp_rmq.query(left, pos + 1, absorb_length)) {
      std::uint64_t low = 0;
      std::uint64_t high = recent_phrase_ends.size();
      std::uint64_t key = text_window_offset +
        (std::uint64_t)sa_of_rev_text_window[left - 1];

      while (high - low > 1) {
        std::uint64_t mid = (low + high) / 2;
        if (text_window_offset + (std::uint64_t)recent_phrase_ends[mid] > key)
          high = mid;
        else low = mid;
      }

      prev_phrase_id = nonrecent_phrase_count + low;
      return true;
    }
  }

  std::uint64_t right = marked.succ(pos);
  if (right < window_size && right == isa_pos_to_ignore)
    right = marked.succ(isa_pos_to_ignore);
  if (right < window_size) {
    if (lcp_rmq.query(pos + 1, right + 1, absorb_length)) {
      std::uint64_t low = 0;
      std::uint64_t high = recent_phrase_ends.size();
      std::uint64_t key = text_window_offset +
        (std::uint64_t)sa_of_rev_text_window[right];

      while (high - low > 1) {
        std::uint64_t mid = (low + high) / 2;
        if (text_window_offset + (std::uint64_t)recent_phrase_ends[mid] > key)
          high = mid;
        else low = mid;
      }

      prev_phrase_id = nonrecent_phrase_count + low;
      return true;
    }
  }

  return false;
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type>
bool absorb_two(
    std::uint64_t text_length,
    std::uint64_t prev_phrase_id,
    std::uint64_t parsed_prefix_length,
    std::uint64_t hash_window_offset,
    std::uint64_t mersenne_prime_exponent,
    const space_efficient_vector<phrase<char_type,
        text_offset_type, phrase_offset_type> > &parsing,
    const space_efficient_vector<std::uint64_t> &phrase_hashes,
    const space_efficient_vector<node_id_type> &leaf_pointers,
    const space_efficient_vector<trie_node<char_type,
        text_offset_type, node_id_type> > &trie_nodes,
    const phrase_offset_type *last_phrase_length,
    const text_offset_type *last_phrase_link_if_in_trie,
    const std::uint64_t *hash_power,
    const std::uint64_t *hash_prefix_of_rev_text,
    std::uint64_t phrase_info_offset) {

  std::uint64_t cur_parsing_size = parsing.size();
  std::uint64_t two_last_phrases_len =
    (std::uint64_t)parsing[cur_parsing_size - 2].m_len +
    (std::uint64_t)parsing[cur_parsing_size - 1].m_len;

  return common_part(text_length, prev_phrase_id, parsed_prefix_length,
      two_last_phrases_len, hash_window_offset, mersenne_prime_exponent,
      parsing, phrase_hashes, leaf_pointers, trie_nodes, last_phrase_length,
      last_phrase_link_if_in_trie, hash_power, hash_prefix_of_rev_text,
      phrase_info_offset);
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type>
bool absorb_one(
    std::uint64_t text_length,
    std::uint64_t prev_phrase_id,
    std::uint64_t parsed_prefix_length,
    std::uint64_t hash_window_offset,
    std::uint64_t mersenne_prime_exponent,
    const space_efficient_vector<phrase<char_type,
        text_offset_type, phrase_offset_type> > &parsing,
    const space_efficient_vector<std::uint64_t> &phrase_hashes,
    const space_efficient_vector<node_id_type> &leaf_pointers,
    const space_efficient_vector<trie_node<char_type,
        text_offset_type, node_id_type> > &trie_nodes,
    const phrase_offset_type *last_phrase_length,
    const text_offset_type *last_phrase_link_if_in_trie,
    const std::uint64_t *hash_power,
    const std::uint64_t *hash_prefix_of_rev_text,
    std::uint64_t phrase_info_offset) {

  std::uint64_t cur_parsing_size = parsing.size();
  std::uint64_t last_phrase_len = parsing[cur_parsing_size - 1].m_len;

  if (parsing[prev_phrase_id].m_len < last_phrase_len)
    return common_part(text_length, prev_phrase_id, parsed_prefix_length,
        last_phrase_len, hash_window_offset, mersenne_prime_exponent,
        parsing, phrase_hashes, leaf_pointers, trie_nodes, last_phrase_length,
        last_phrase_link_if_in_trie, hash_power, hash_prefix_of_rev_text,
        phrase_info_offset);

  if (parsing[prev_phrase_id].m_char != parsing.back().m_char)
    return false;
  else if (last_phrase_len == 1)
    return true;
  else if (last_phrase_link_if_in_trie[
      parsed_prefix_length - 1 - phrase_info_offset] ==
      std::numeric_limits<text_offset_type>::max())
    return false;

  return (std::uint64_t)trie_nodes[nearest_common_ancestor(
      (std::uint64_t)
        leaf_pointers[(std::uint64_t)last_phrase_link_if_in_trie[
        parsed_prefix_length - 1 - phrase_info_offset]],
      (std::uint64_t)
        leaf_pointers[(std::uint64_t)parsing[prev_phrase_id].m_link],
        trie_nodes)].m_depth + 1 >= last_phrase_len;
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type>
bool common_part(
    std::uint64_t text_length,
    std::uint64_t prev_phrase_id,
    std::uint64_t parsed_prefix_length,
    std::uint64_t substring_length,
    std::uint64_t hash_window_offset,
    std::uint64_t mersenne_prime_exponent,
    const space_efficient_vector<phrase<
        char_type, text_offset_type, phrase_offset_type> > &parsing,
    const space_efficient_vector<std::uint64_t> &phrase_hashes,
    const space_efficient_vector<node_id_type> &leaf_pointers,
    const space_efficient_vector<trie_node<char_type,
        text_offset_type, node_id_type> > &trie_nodes,
    const phrase_offset_type *last_phrase_length,
    const text_offset_type *last_phrase_link_if_in_trie,
    const std::uint64_t *hash_power,
    const std::uint64_t *hash_prefix_of_rev_text,
    std::uint64_t phrase_info_offset) {

  if ((std::uint64_t)parsing[prev_phrase_id].m_len >= substring_length ||
      phrase_hashes[prev_phrase_id] !=
      hash_of_reversed_substr<char_type>(parsed_prefix_length -
        (std::uint64_t)parsing[prev_phrase_id].m_len,
        parsed_prefix_length, text_length, hash_window_offset,
        mersenne_prime_exponent, hash_power, hash_prefix_of_rev_text))
    return false;

  std::uint64_t pos = parsed_prefix_length - parsing[prev_phrase_id].m_len;
  std::uint64_t pos_offset = pos - phrase_info_offset;
  if ((std::uint64_t)last_phrase_length[pos_offset] +
      (std::uint64_t)parsing[prev_phrase_id].m_len !=
      substring_length + 1 ||
      last_phrase_link_if_in_trie[pos_offset] ==
      std::numeric_limits<text_offset_type>::max() ||
      prev_phrase_id == 0)
    return false;

  return (std::uint64_t)trie_nodes[nearest_common_ancestor(
      (std::uint64_t)
        leaf_pointers[(std::uint64_t)last_phrase_link_if_in_trie[pos_offset]],
      (std::uint64_t)
        leaf_pointers[prev_phrase_id - 1], trie_nodes)].m_depth +
        (std::uint64_t)parsing[prev_phrase_id].m_len >= substring_length;
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type,
  typename window_offset_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose) {

  utils::initialize_stats();
  srand(time(0) + getpid());
  std::uint64_t text_length = utils::file_size(text_filename);

  if (text_length == 0) {
    fprintf(stderr, "\nError: text length is 0.\n");
    std::exit(EXIT_FAILURE);
  }

  if (phrase_length_limit == 0) {
    fprintf(stderr, "\nError: phrase length limit is 0.\n");
    std::exit(EXIT_FAILURE);
  }

  if ((std::uint64_t)std::numeric_limits<text_offset_type>::max() <
      text_length - 1) {
    fprintf(stderr, "\nError: sizeof(text_offset_type) = %lu "
        "is too small.\n", sizeof(text_offset_type));
    std::exit(EXIT_FAILURE);
  }

  if ((std::uint64_t)std::numeric_limits<node_id_type>::max() <
      2 * text_length - 1) {
    fprintf(stderr, "\nError: sizeof(node_id_type) = %lu "
        "is too small.\n", sizeof(node_id_type));
    std::exit(EXIT_FAILURE);
  }

  if ((std::uint64_t)std::numeric_limits<phrase_offset_type>::max() <
      phrase_length_limit - 1) {
    fprintf(stderr, "\nError: sizeof(phrase_offset_type) = %lu "
        "is too small.\n", sizeof(phrase_offset_type));
    std::exit(EXIT_FAILURE);
  }

  if ((std::uint64_t)std::numeric_limits<window_offset_type>::max() <
      3 * phrase_length_limit - 1) {
    fprintf(stderr, "\nError: sizeof(window_offset_type) = %lu "
        "is too small.\n", sizeof(window_offset_type));
    std::exit(EXIT_FAILURE);
  }

  typedef trie_node<char_type, text_offset_type, node_id_type> node_type;

  long double parsing_start = utils::wclock();

  std::FILE *file_text =
    utils::file_open(text_filename, "r");

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "Running LZ-End Toolkit v0.1.0\n");
  fprintf(stderr, "Mode = parsing\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length,
      (1.L * text_length * sizeof(char_type)) / (1UL << 20));
  fprintf(stderr, "Phrase length limit = %lu\n", phrase_length_limit);
  fprintf(stderr, "sizeof(char_type) = %lu\n",
      sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "Verbose = %s\n", verbose ? "TRUE" : "FALSE");
  fprintf(stderr, "\n\n");

  // Initialize the arrays used to compute hashes of reversed substrings.
  static const std::uint64_t mersenne_prime_exponent = 61;
  std::uint64_t hash_variable =
    rand_mod_mersenne(mersenne_prime_exponent);

  std::uint64_t max_needed_hash_length =
    std::min(text_length, phrase_length_limit * 2);
  std::uint64_t *hash_power =
    utils::allocate_array<std::uint64_t>(max_needed_hash_length + 1);
  std::uint64_t *hash_prefix_of_rev_text_window =
    utils::allocate_array<std::uint64_t>(max_needed_hash_length + 1);
  hash_power[0] = 1;
  for (std::uint64_t i = 1; i <= max_needed_hash_length; ++i)
    hash_power[i] = mul_mod_meresenne(hash_power[i - 1],
        hash_variable, mersenne_prime_exponent);

  // Allocate the arrays storing info about the
  // last phrase in the parsing of every prefix.
  phrase_offset_type *last_phrase_length =
    utils::allocate_array<phrase_offset_type>(max_needed_hash_length);
  text_offset_type *last_phrase_link_if_in_trie =
    utils::allocate_array<text_offset_type>(max_needed_hash_length);

  // Allocate text window.
  std::uint64_t max_text_window_size =
    std::min(text_length, phrase_length_limit * 3);
  char_type *text_window =
    utils::allocate_array<char_type>(max_text_window_size);

  // Allocate the root of the trie.
  typedef space_efficient_vector<node_type> trie_nodes_vector_type;
  trie_nodes_vector_type *trie_nodes = new trie_nodes_vector_type();
  trie_nodes->push_back(node_type());
  std::uint64_t trie_root = 0;
  (*trie_nodes)[trie_root].m_depth = 0;

  // Create the hash_table used to navigate the trie.
  typedef packed_pair<phrase_offset_type, std::uint64_t> nav_key_type;
  typedef node_id_type nav_value_type;
  typedef hash_table<nav_key_type, nav_value_type, node_id_type> nav_type;
  nav_type *nav = new nav_type();

  // Allocate suffix array, inverse suffix array,
  // LCP array and predecessor tree.
  window_offset_type *sa_of_rev_text_window =
    utils::allocate_array<window_offset_type>(max_text_window_size);
  window_offset_type *isa_of_rev_text_window =
    utils::allocate_array<window_offset_type>(max_text_window_size);
  window_offset_type *lcp_array =
    utils::allocate_array<window_offset_type>(max_text_window_size);
  predecessor_tree *marked = new predecessor_tree(max_text_window_size);

  // Process blocks left-to-right.
  typedef space_efficient_vector<node_id_type> leaf_pointers_vector_type;
  typedef phrase<char_type, text_offset_type, phrase_offset_type> phrase_type;
  typedef space_efficient_vector<phrase_type> parsing_vector_type;
  leaf_pointers_vector_type *leaf_pointers = new leaf_pointers_vector_type();
  parsing_vector_type *parsing = new parsing_vector_type();
  space_efficient_vector<std::uint64_t> *phrase_hashes =
    new space_efficient_vector<std::uint64_t>();
  std::uint64_t n_blocks =
    (text_length + phrase_length_limit - 1) / phrase_length_limit;
  std::uint64_t longest_prefix_in_trie = 0;
  std::uint64_t max_permanent_phrase_length = 0;
  std::uint64_t block_group_size = (n_blocks + 9999) / 10000;
  for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {

    // Compute basic parameters.
    std::uint64_t block_beg = block_id * phrase_length_limit;
    std::uint64_t block_end =
      std::min(block_beg + phrase_length_limit, text_length);
    std::uint64_t block_size = block_end - block_beg;

    // Print info about the block.
    if (verbose) {
      fprintf(stderr, "Process block %lu/%lu [%lu..%lu):\n",
          block_id + 1, n_blocks, block_beg, block_end);
    } else if (block_id % block_group_size == 0)
      fprintf(stderr, "\rParse: %.2Lf%%",
          ((100.L * block_id + 1)) / n_blocks);

    // Prepare window of text.
    std::uint64_t text_window_offset = 0;
    if (block_id > 2) {
      text_window_offset = (block_id - 2) * phrase_length_limit;

      // Shift the symbols from the previous window to the left.
      for (std::uint64_t i = 0; i < 2 * phrase_length_limit; ++i)
        text_window[i] = text_window[phrase_length_limit + i];
    }
    std::uint64_t text_window_size = block_end - text_window_offset;

    // Read the current block into the window.
    {
      if (verbose)
        fprintf(stderr, "  Read block: ");
      long double read_start = utils::wclock();

      utils::read_from_file(text_window + block_beg -
          text_window_offset, block_size, file_text);

      long double elapsed = utils::wclock() - read_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    // Compute Karp-Rabin hashes.
    std::uint64_t hash_window_offset = text_length - block_end;
    std::uint64_t hash_window_size =
      std::min(phrase_length_limit * 2, block_end) + 1;
    {
      if (verbose)
        fprintf(stderr, "  Compute hashes: ");
      long double hash_start = utils::wclock();

      hash_prefix_of_rev_text_window[0] = 0;
      for (std::uint64_t i = 1; i < hash_window_size; ++i) {
        std::uint64_t ab_mod_p = mul_mod_meresenne(
            hash_prefix_of_rev_text_window[i - 1],
            hash_variable, mersenne_prime_exponent);
        std::uint64_t c_mod_p = mod_mersenne(
            (std::uint64_t)text_window[block_end - i - text_window_offset],
            mersenne_prime_exponent);
        hash_prefix_of_rev_text_window[i] = mod_mersenne(
            ab_mod_p + c_mod_p, mersenne_prime_exponent);
      }

      long double elapsed = utils::wclock() - hash_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    // Prepare window with info about recent phrases.
    std::uint64_t phrase_info_offset = 0;
    if (block_id > 1) {
      phrase_info_offset = (block_id - 1) * phrase_length_limit;

      // Shift the data from previous window to the left.
      for (std::uint64_t i = 0; i < phrase_length_limit; ++i) {
        last_phrase_length[i] = last_phrase_length[i + phrase_length_limit];
        last_phrase_link_if_in_trie[i] =
          last_phrase_link_if_in_trie[i + phrase_length_limit];
      }
    }

    // Prepare SA, ISA, and LCP of current text window reversed.
    // Note: sais uses n integers of extra space, but since our
    // LCP array construction also allocated n extra integers, sais
    // does not increase the peak RAM allocation of the algorithm.
    {
      if (verbose)
        fprintf(stderr, "  Compute SA: ");
      long double compute_sa_start = utils::wclock();

      std::reverse(text_window, text_window + text_window_size);
      construct_sa(text_window, text_window_size, sa_of_rev_text_window);

      long double elapsed = utils::wclock() - compute_sa_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    {
      if (verbose)
        fprintf(stderr, "  Compute LCP: ");
      long double compute_lcp_start = utils::wclock();

      construct_lcp(text_window, text_window_size,
          sa_of_rev_text_window, lcp_array);
      std::reverse(text_window, text_window + text_window_size);

      long double elapsed = utils::wclock() - compute_lcp_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    {
      if (verbose)
        fprintf(stderr, "  Compute ISA: ");
      long double compute_isa_start = utils::wclock();

      for (std::uint64_t i = 0; i < text_window_size; ++i)
        sa_of_rev_text_window[i] = text_window_size - 1 -
          (std::uint64_t)sa_of_rev_text_window[i];
      for (std::uint64_t i = 0; i < text_window_size; ++i) {
        std::uint64_t addr = sa_of_rev_text_window[i];
        isa_of_rev_text_window[addr] = i;
      }

      long double elapsed = utils::wclock() - compute_isa_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    // Create the list of recent phrase ends
    // and prepare the marked array.
    space_efficient_vector<window_offset_type> recent_phrase_ends;
    {
      if (verbose)
        fprintf(stderr, "  Mark recent: ");
      long double mark_start = utils::wclock();

      marked->clear();
      std::uint64_t total_len = longest_prefix_in_trie;
      for (std::uint64_t phrase_id = leaf_pointers->size();
          phrase_id < parsing->size(); ++phrase_id) {
        std::uint64_t length = (*parsing)[phrase_id].m_len;
        std::uint64_t end = total_len + length - 1;
        std::uint64_t offset = end - text_window_offset;
        std::uint64_t sa_idx = isa_of_rev_text_window[offset];
        if (phrase_id + 1 < parsing->size())
          marked->set(sa_idx);
        recent_phrase_ends.push_back((window_offset_type)offset);
        total_len += length;
      }

      long double elapsed = utils::wclock() - mark_start;
      if (verbose)
        fprintf(stderr, "%.2Lfs\n", elapsed);
    }

    // Parse current block.
    parse_block(text_length, phrase_length_limit, block_id,
        mersenne_prime_exponent, hash_variable, hash_window_offset,
        longest_prefix_in_trie, *parsing, *phrase_hashes,
        *leaf_pointers, recent_phrase_ends, trie_root, *trie_nodes,
        last_phrase_length, last_phrase_link_if_in_trie, hash_power,
        hash_prefix_of_rev_text_window, *nav, text_window, text_window_offset,
        phrase_info_offset, sa_of_rev_text_window, isa_of_rev_text_window,
        lcp_array, verbose, *marked, max_permanent_phrase_length);

    long double elapsed = utils::wclock() - parsing_start;
    if (verbose) {
      fprintf(stderr, "  Cur parsing size = %lu\n", parsing->size());
      fprintf(stderr, "  Cur avg phrase length = %.2Lf\n",
          (1.L * block_end) / parsing->size());
      fprintf(stderr, "  Max permanent phrase length = %lu\n",
          max_permanent_phrase_length);
      fprintf(stderr, "  Peak RAM usage so far = %.2LfMiB\n",
          (1.L * utils::get_peak_ram_allocation()) / (1L << 20));
      fprintf(stderr, "  Total time so far = %.2Lfs (%.2Lfus/symbol)\n",
          elapsed, (elapsed * 1000000.L) / block_end);
    }
  }

  if (!verbose)
    fprintf(stderr, "\rParse: 100.00%%\n");

  // Clean up.
  std::fclose(file_text);
  delete marked;
  utils::deallocate(lcp_array);
  utils::deallocate(isa_of_rev_text_window);
  utils::deallocate(sa_of_rev_text_window);
  delete nav;
  delete trie_nodes;
  utils::deallocate(text_window);
  utils::deallocate(last_phrase_link_if_in_trie);
  utils::deallocate(last_phrase_length);
  utils::deallocate(hash_prefix_of_rev_text_window);
  utils::deallocate(hash_power);
  delete phrase_hashes;
  delete leaf_pointers;

  // Determine basic statistics.
  std::uint64_t parsing_size = parsing->size();
  std::uint64_t max_phrase_len = 0;
  for (std::uint64_t i = 0; i < parsing->size(); ++i)
    max_phrase_len =
      std::max(max_phrase_len, (std::uint64_t)(*parsing)[i].m_len);

  // Write parsing to file.
  {

    // Allocate buffer.
    static const std::uint64_t bufsize = (64 << 10);
    typedef phrase<char_type, text_offset_type, text_offset_type>
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

  // Clean up.
  delete parsing;

  // Print summary.
  long double elapsed = utils::wclock() - parsing_start;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  Absolute time = %.2Lfs\n", elapsed);
  fprintf(stderr, "  Relative time = %.2Lfus/symbol\n",
      (elapsed * 1000000.L) / text_length);
  fprintf(stderr, "  Number of phrases = %lu\n", parsing_size);
  fprintf(stderr, "  Maximal phrase length = %lu\n", max_phrase_len);
  fprintf(stderr, "  Average phrase length = %.2Lf\n",
      (1.L * text_length) / parsing_size);
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type,
  typename window_offset_type>
void parse_block(
    std::uint64_t text_length,
    std::uint64_t phrase_length_limit,
    std::uint64_t block_id,
    std::uint64_t mersenne_prime_exponent,
    std::uint64_t hash_variable,
    std::uint64_t hash_window_offset,
    std::uint64_t &longest_prefix_in_trie,
    space_efficient_vector<phrase<char_type,
        text_offset_type, phrase_offset_type> > &parsing,
    space_efficient_vector<std::uint64_t> &phrase_hashes,
    space_efficient_vector<node_id_type> &leaf_pointers,
    space_efficient_vector<window_offset_type> &recent_phrase_ends,
    std::uint64_t trie_root,
    space_efficient_vector<trie_node<
        char_type,text_offset_type, node_id_type> > &trie_nodes,
    phrase_offset_type *last_phrase_length,
    text_offset_type *last_phrase_link_if_in_trie,
    std::uint64_t *hash_power,
    std::uint64_t *hash_prefix_of_rev_text,
    hash_table<packed_pair<phrase_offset_type,
        std::uint64_t>, node_id_type, node_id_type> &nav,
    const char_type *text_window,
    std::uint64_t text_window_offset,
    std::uint64_t phrase_info_offset,
    const window_offset_type *sa_of_rev_text_window,
    const window_offset_type *isa_of_rev_text_window,
    const window_offset_type *lcp_array,
    bool verbose,
    predecessor_tree &marked,
    std::uint64_t &max_permanent_phrase_length) {

  std::uint64_t n_blocks =
    (text_length + phrase_length_limit - 1) / phrase_length_limit;
  std::uint64_t block_beg = block_id * phrase_length_limit;
  std::uint64_t block_end =
    std::min(block_beg + phrase_length_limit, text_length);
  std::uint64_t window_size = block_end - text_window_offset;
  std::uint64_t nonrecent_phrase_count = leaf_pointers.size();

  long double start = utils::wclock();
  if (verbose)
    fprintf(stderr, "  Parse block:\n");

  typedef trie_node<char_type, text_offset_type, node_id_type> node_type;
  typedef phrase<char_type, text_offset_type, phrase_offset_type> phrase_type;
  typedef packed_pair<phrase_offset_type, std::uint64_t> nav_key_type;
  typedef node_id_type nav_value_type;
  typedef rmq_tree<window_offset_type, 64> rmq_type;


  rmq_type *lcp_rmq = NULL;
  {
    if (verbose)
      fprintf(stderr, "    Initialize RMQ: ");
    long double rmq_init_start = utils::wclock();
    lcp_rmq = new rmq_type(lcp_array, window_size);
    long double elapsed = utils::wclock() - rmq_init_start;
    if (verbose)
      fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  if (verbose)
    fprintf(stderr, "    Parse: ");
  long double parsing_start = utils::wclock();

  for (std::uint64_t parsed_prefix_length = block_beg;
      parsed_prefix_length < block_end; ++parsed_prefix_length) {

    // Extend the parsing to include text[i].
    std::uint64_t cur_parsing_size = parsing.size();
    std::uint64_t prev_phrase_id = 0;
    last_phrase_link_if_in_trie[parsed_prefix_length - phrase_info_offset]
      = std::numeric_limits<text_offset_type>::max();

    if (cur_parsing_size >= 1) {

      // Compute the phrase candidate.
      {
        if (cur_parsing_size >= 2 &&
            (std::uint64_t)parsing[cur_parsing_size - 2].m_len +
            (std::uint64_t)parsing[cur_parsing_size - 1].m_len + 1
              <= phrase_length_limit) {

          // Note that the candidate found below is useful even
          // when checking if the last phrase had a previous occurrence.
          std::uint64_t two_last_phrases_len =
            parsing[cur_parsing_size - 2].m_len +
            parsing[cur_parsing_size - 1].m_len;
          prev_phrase_id = approx_find(text_window, text_window_offset,
              parsed_prefix_length, two_last_phrases_len, trie_root,
              text_length, hash_window_offset, mersenne_prime_exponent,
              hash_power, hash_prefix_of_rev_text, trie_nodes, nav);
        } else if ((std::uint64_t)parsing[cur_parsing_size - 1].m_len + 1
            <= phrase_length_limit) {

          std::uint64_t last_phrase_len = parsing[cur_parsing_size - 1].m_len;
          prev_phrase_id = approx_find(text_window, text_window_offset,
              parsed_prefix_length, last_phrase_len, trie_root,
              text_length, hash_window_offset, mersenne_prime_exponent,
              hash_power, hash_prefix_of_rev_text, trie_nodes, nav);
        }
      }

      // Try absorbing the last two phrases.
      bool found = false;
      if (cur_parsing_size >= 2) {
        std::uint64_t last_phrase_len = parsing[cur_parsing_size - 1].m_len;
        std::uint64_t two_last_phrases_len = last_phrase_len +
          (std::uint64_t)parsing[cur_parsing_size - 2].m_len;

        if (trie_nodes[trie_root].has_child() &&
            two_last_phrases_len + 1 <= phrase_length_limit &&
            absorb_two<char_type, text_offset_type, node_id_type,
              phrase_offset_type>(text_length, prev_phrase_id,
                parsed_prefix_length, hash_window_offset,
                mersenne_prime_exponent, parsing, phrase_hashes,
                leaf_pointers, trie_nodes, last_phrase_length,
                last_phrase_link_if_in_trie, hash_power,
                hash_prefix_of_rev_text, phrase_info_offset)) {
          std::uint64_t text_pos_to_reset =
            parsed_prefix_length - last_phrase_len - 1;
          std::uint64_t isa_pos_to_reset =
            isa_of_rev_text_window[text_pos_to_reset - text_window_offset];
          marked.reset(isa_pos_to_reset);
          parsing.pop_back();
          phrase_hashes.pop_back();
          recent_phrase_ends.pop_back();
          parsing.back().m_len = two_last_phrases_len + 1;
          last_phrase_link_if_in_trie[
            parsed_prefix_length - phrase_info_offset] = prev_phrase_id;
          found = true;
        } else if (two_last_phrases_len + 1 <= phrase_length_limit) {
          std::uint64_t text_pos_to_ignore =
            parsed_prefix_length - last_phrase_len - 1;
          std::uint64_t isa_pos_to_ignore =
            isa_of_rev_text_window[text_pos_to_ignore - text_window_offset];
          if (marked_lcp<char_type, text_offset_type, phrase_offset_type,
              window_offset_type>(parsed_prefix_length - 1,
                text_window_offset, window_size, two_last_phrases_len,
                isa_pos_to_ignore, sa_of_rev_text_window,
                isa_of_rev_text_window, marked, recent_phrase_ends,
                nonrecent_phrase_count, prev_phrase_id, *lcp_rmq)) {
            marked.reset(isa_pos_to_ignore);
            parsing.pop_back();
            phrase_hashes.pop_back();
            recent_phrase_ends.pop_back();
            parsing.back().m_len = two_last_phrases_len + 1;
            found = true;
          }
        }
      }

      // Try absorbing the last phrase.
      if (!found) {
        std::uint64_t last_phrase_len = parsing[cur_parsing_size - 1].m_len;
        if (trie_nodes[trie_root].has_child() &&
            last_phrase_len + 1 <= phrase_length_limit &&
            absorb_one(text_length, prev_phrase_id, parsed_prefix_length,
              hash_window_offset, mersenne_prime_exponent, parsing,
              phrase_hashes, leaf_pointers, trie_nodes, last_phrase_length,
              last_phrase_link_if_in_trie, hash_power, hash_prefix_of_rev_text,
              phrase_info_offset)) {
          parsing.back().m_len = last_phrase_len + 1;
          last_phrase_link_if_in_trie[
            parsed_prefix_length - phrase_info_offset] = prev_phrase_id;
          found = true;
        } else if (last_phrase_len + 1 <= phrase_length_limit) {
          if (marked_lcp<char_type, text_offset_type, phrase_offset_type,
              window_offset_type>(parsed_prefix_length - 1,
                text_window_offset, window_size, last_phrase_len,
                window_size, sa_of_rev_text_window, isa_of_rev_text_window,
                marked, recent_phrase_ends, nonrecent_phrase_count,
                prev_phrase_id, *lcp_rmq)) {
            parsing.back().m_len = last_phrase_len + 1;
            found = true;
          }
        }
      }

      // Do not absorb any phrases.
      if (!found) {
        parsing.push_back(phrase_type());
        parsing.back().m_len = 1;
        phrase_hashes.push_back(std::uint64_t());
        recent_phrase_ends.push_back(window_offset_type());
      }
    } else {
      parsing.push_back(phrase_type());
      parsing.back().m_len = 1;
      phrase_hashes.push_back(std::uint64_t());
      recent_phrase_ends.push_back(window_offset_type());
    }


    if ((std::uint64_t)parsing.back().m_len == 1 && parsing.size() > 1) {
      std::uint64_t text_pos_to_set =
        parsed_prefix_length - 1;
      std::uint64_t isa_pos_to_set =
        isa_of_rev_text_window[text_pos_to_set - text_window_offset];
      marked.set(isa_pos_to_set);
    }
    last_phrase_length[parsed_prefix_length - phrase_info_offset] =
      parsing.back().m_len;
    parsing.back().m_char =
      text_window[parsed_prefix_length - text_window_offset];
    phrase_hashes.back() = hash_of_reversed_substr<char_type>(
        parsed_prefix_length - parsing.back().m_len + 1,
        parsed_prefix_length + 1, text_length, hash_window_offset,
        mersenne_prime_exponent, hash_power, hash_prefix_of_rev_text);
    parsing.back().m_link = prev_phrase_id;
    recent_phrase_ends.back() = parsed_prefix_length - text_window_offset;
  }

  long double parse_time = utils::wclock() - parsing_start;
  if (verbose)
    fprintf(stderr, "%.2Lfs\n", parse_time);

  if (block_id + 1 < n_blocks) {
    if (verbose)
      fprintf(stderr, "    Update trie: ");
    long double update_trie_start = utils::wclock();

    // Update the trie.
    std::uint64_t old_leaf_pointers_size = leaf_pointers.size();
    std::uint64_t old_longest_prefix_in_trie = longest_prefix_in_trie;
    std::uint64_t local_max_permanent_phrase_length = 0;
    for (std::uint64_t phrase_id = old_leaf_pointers_size;
        phrase_id < parsing.size() && longest_prefix_in_trie <= block_beg;
        ++phrase_id) {
      longest_prefix_in_trie += (std::uint64_t)parsing[phrase_id].m_len;
      local_max_permanent_phrase_length =
        std::max(local_max_permanent_phrase_length,
            (std::uint64_t)parsing[phrase_id].m_len);

      // Descend 'blindly' in the tree.
      std::uint64_t cur_node = trie_root;
      while ((std::uint64_t)trie_nodes[cur_node].m_depth <
          (std::uint64_t)parsing[phrase_id].m_len) {
        std::uint64_t cur_depth = trie_nodes[cur_node].m_depth;
        std::uint64_t next_char_idx = longest_prefix_in_trie - 1 - cur_depth;
        char_type next_char = text_window[next_char_idx - text_window_offset];
        std::uint64_t child =
          trie_nodes[cur_node].get_child(next_char, trie_nodes);
        if (child == std::numeric_limits<std::uint64_t>::max()) break;
        else cur_node = child;
      }

      // Go up in the tree to find the
      // insertion node or it's parent.
      std::uint64_t lcs_length = 0;
      char_type mismatch_char = 0;
      if (trie_nodes[trie_root].has_child()) {
        lcs_length = longest_common_suffix(text_window, text_window_offset,
            longest_prefix_in_trie, phrase_length_limit,
            (std::uint64_t)trie_nodes[cur_node].m_phr,
            parsing, mismatch_char);
        while (cur_node != trie_root &&
            (std::uint64_t)trie_nodes[trie_nodes[cur_node].m_parent].m_depth
            >= lcs_length) cur_node = trie_nodes[cur_node].m_parent;
      }

      if (lcs_length >= phrase_length_limit) {

        // Special case, we don't insert this prefix into the trie,
        // because it would not affect the relevant part of the trie
        // Instead we redirect left pointer to the other leaf.
        std::uint64_t node_ptr = trie_nodes[cur_node].m_phr;
        leaf_pointers.push_back(leaf_pointers[node_ptr]);
        continue;
      }

      // Insert the new leaf and possibly one new internal node.
      if ((std::uint64_t)trie_nodes[cur_node].m_depth == lcs_length) {
        std::uint64_t new_leaf = trie_nodes.size();
        trie_nodes.push_back(node_type());
        trie_nodes[new_leaf].m_depth = longest_prefix_in_trie;
        trie_nodes[new_leaf].m_phr = phrase_id;
        trie_nodes[new_leaf].m_char = text_window[
          longest_prefix_in_trie - 1 - lcs_length - text_window_offset];
        trie_nodes[new_leaf].m_parent = cur_node;
        trie_nodes[cur_node].add_child(new_leaf, trie_nodes);
        trie_nodes[cur_node].m_phr = phrase_id;  // may be undefined
        leaf_pointers.push_back((node_id_type)new_leaf);

        // Update nav.
        if (cur_node != trie_root &&
            trie_nodes[cur_node].has_one_child(trie_nodes)) {

          std::uint64_t par_node = trie_nodes[cur_node].m_parent;
          std::uint64_t pv = compute_pv(
              (std::uint64_t)trie_nodes[cur_node].m_depth,
              (std::uint64_t)trie_nodes[par_node].m_depth);

          std::uint64_t hv =
            hash_of_reversed_substr<char_type>(
              longest_prefix_in_trie - pv,
              longest_prefix_in_trie,
              text_length, hash_window_offset,
              mersenne_prime_exponent,
              hash_power, hash_prefix_of_rev_text);
          nav.insert(
              nav_key_type((phrase_offset_type)pv, hv),
              nav_value_type(cur_node));
        }
      } else {

        // Split an edge.
        std::uint64_t parent = trie_nodes[cur_node].m_parent;
        std::uint64_t new_internal_node = trie_nodes.size();
        trie_nodes.push_back(node_type());
        trie_nodes[new_internal_node].m_depth = lcs_length;
        trie_nodes[new_internal_node].m_phr = phrase_id;
        trie_nodes[new_internal_node].m_parent = parent;
        trie_nodes[new_internal_node].m_char = trie_nodes[cur_node].m_char;
        trie_nodes[new_internal_node].add_child(cur_node, trie_nodes);
        trie_nodes[parent].remove_child(cur_node, trie_nodes);
        trie_nodes[parent].add_child(new_internal_node, trie_nodes);
        trie_nodes[cur_node].m_parent = new_internal_node;
        trie_nodes[cur_node].m_char = mismatch_char;
        std::uint64_t new_leaf = trie_nodes.size();
        trie_nodes.push_back(node_type());
        trie_nodes[new_leaf].m_depth = longest_prefix_in_trie;
        trie_nodes[new_leaf].m_phr = phrase_id;
        trie_nodes[new_leaf].m_char = text_window[
          longest_prefix_in_trie - 1 - lcs_length - text_window_offset];
        trie_nodes[new_leaf].m_parent = new_internal_node;
        leaf_pointers.push_back((node_id_type)new_leaf);
        trie_nodes[new_internal_node].add_child(new_leaf, trie_nodes);

        // Update nav.
        {
          std::uint64_t pv = compute_pv(
              (std::uint64_t)trie_nodes[new_internal_node].m_depth,
              (std::uint64_t)trie_nodes[parent].m_depth);
          std::uint64_t hv =
            hash_of_reversed_substr<char_type>(
              longest_prefix_in_trie - pv,
              longest_prefix_in_trie,
              text_length, hash_window_offset,
              mersenne_prime_exponent,
              hash_power, hash_prefix_of_rev_text);

          nav_key_type key((phrase_offset_type)pv, hv);
          nav_value_type value = new_internal_node;
          nav_value_type *ret = nav.find(key);
          if (ret == NULL) nav.insert(key, value);
          else *ret = new_internal_node;
        }

        if (trie_nodes[cur_node].has_child()) {
          std::uint64_t oldpv = compute_pv(
              (std::uint64_t)trie_nodes[cur_node].m_depth,
              (std::uint64_t)trie_nodes[parent].m_depth);
          if (oldpv <= (std::uint64_t)trie_nodes[new_internal_node].m_depth) {
            std::uint64_t pv = compute_pv(
                (std::uint64_t)trie_nodes[cur_node].m_depth,
                (std::uint64_t)trie_nodes[new_internal_node].m_depth);
            std::uint64_t hv = hash_of_reversed_substr(
                (std::uint64_t)trie_nodes[cur_node].m_phr,
                pv, mersenne_prime_exponent, hash_variable, parsing);
            nav.insert(
                nav_key_type((phrase_offset_type)pv, hv),
                nav_value_type(cur_node));
          }
        }
      }
    }

    max_permanent_phrase_length =
      std::max(max_permanent_phrase_length,
          local_max_permanent_phrase_length);

    if (verbose) {
      long double update_trie_time = utils::wclock() - update_trie_start;
      fprintf(stderr, "%.2Lfs\n", update_trie_time);
    }

    if (verbose)
      fprintf(stderr, "    Update links: ");
    long double phrase_pointer_update_start = utils::wclock();

    // Compute marked positions for strings newly inserted into trie,
    // except the last one, which overlap the current block.
    marked.clear();
    recent_phrase_ends.clear();
    std::uint64_t new_leaf_pointers_size = leaf_pointers.size();
    std::uint64_t total_len = old_longest_prefix_in_trie;
    for (std::uint64_t i = old_leaf_pointers_size;
        i + 1 < new_leaf_pointers_size; ++i) {
      std::uint64_t length = parsing[i].m_len;
      std::uint64_t end = total_len + length - 1;
      std::uint64_t offset = end - text_window_offset;
      std::uint64_t sa_idx = isa_of_rev_text_window[offset];
      marked.set(sa_idx);
      recent_phrase_ends.push_back((window_offset_type)offset);
      total_len += length;
    }

    // Update last_phrase_link_if_in_trie for the current block.
    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      if ((std::uint64_t)last_phrase_length[i - phrase_info_offset] > 1 &&
          last_phrase_link_if_in_trie[i - phrase_info_offset] ==
          std::numeric_limits<text_offset_type>::max()) {

        std::uint64_t prev_phrase_id = 0;
        if (marked_lcp<char_type, text_offset_type, phrase_offset_type,
            window_offset_type>(i - 1, text_window_offset, window_size,
              (std::uint64_t)last_phrase_length[i - phrase_info_offset] - 1,
              window_size, sa_of_rev_text_window, isa_of_rev_text_window,
              marked, recent_phrase_ends, nonrecent_phrase_count,
              prev_phrase_id, *lcp_rmq))

          last_phrase_link_if_in_trie[i - phrase_info_offset] =
            prev_phrase_id;
      }

      // Mark the position corresponding to the
      // latest prefix inserted into the trie.
      if (i == total_len +
          (std::uint64_t)parsing[new_leaf_pointers_size - 1].m_len - 1) {
        std::uint64_t length = parsing[new_leaf_pointers_size - 1].m_len;
        std::uint64_t end = total_len + length - 1;
        std::uint64_t offset = end - text_window_offset;
        std::uint64_t sa_idx = isa_of_rev_text_window[offset];
        marked.set(sa_idx);
        recent_phrase_ends.push_back((window_offset_type)offset);
      }
    }

    long double phrase_pointer_update_time = utils::wclock() -
      phrase_pointer_update_start;
    if (verbose)
      fprintf(stderr, "%.2Lfs\n", phrase_pointer_update_time);
  }

  // Clean up.
  delete lcp_rmq;

  // Print summary.
  long double total_time = utils::wclock() - start;
  if (verbose)
    fprintf(stderr, "    Total time = %.2Lfs\n", total_time);
}

template<typename char_type,
  typename text_offset_type,
  typename node_id_type,
  typename phrase_offset_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose) {

  // The condition below could in principle be <= (1 << 32),
  // but sais (suffix sorting algorithm) requires < (1 << 31).
  std::uint64_t phrase_length_limit3 = 3 * phrase_length_limit;
  if (phrase_length_limit3 < ((std::uint64_t)1 << 31))
    parse<char_type, text_offset_type, node_id_type,
      phrase_offset_type, std::uint32_t>(text_filename,
          output_filename, phrase_length_limit, verbose);
  else
    parse<char_type, text_offset_type, node_id_type,
      phrase_offset_type, std::uint64_t>(text_filename,
          output_filename, phrase_length_limit, verbose);
}

template<typename char_type,
  typename text_offset_type,
  typename phrase_offset_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose) {

  std::uint64_t text_length =
    utils::file_size(text_filename) / sizeof(char_type);
  std::uint64_t text_length2 = (text_length << 1);
  if (text_length2 <= ((std::uint64_t)1 << 32))
    parse<char_type, text_offset_type, std::uint32_t, phrase_offset_type>(
        text_filename, output_filename, phrase_length_limit, verbose);
  else if (text_length2 <= ((std::uint64_t)1 << 40))
    parse<char_type, text_offset_type, uint40, phrase_offset_type>(
        text_filename, output_filename, phrase_length_limit, verbose);
  else if (text_length2 <= ((std::uint64_t)1 << 48))
    parse<char_type, text_offset_type, uint48, phrase_offset_type>(
        text_filename, output_filename, phrase_length_limit, verbose);
  else
    parse<char_type, text_offset_type, std::uint64_t, phrase_offset_type>(
        text_filename, output_filename, phrase_length_limit, verbose);
}

template<typename char_type,
  typename text_offset_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose) {

  if (phrase_length_limit <= ((std::uint64_t)1 << 24))
    parse<char_type, text_offset_type, uint24>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else if (phrase_length_limit <= ((std::uint64_t)1 << 32))
    parse<char_type, text_offset_type, std::uint32_t>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else
    parse<char_type, text_offset_type, std::uint64_t>(text_filename,
        output_filename, phrase_length_limit, verbose);
}

}  // namespace lz_end_toolkit_private

template<typename char_type,
  typename text_offset_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose) {
  lz_end_toolkit_private::parse<char_type, text_offset_type>(
      text_filename, output_filename, phrase_length_limit, verbose);
}

#endif  // __LZ_END_TOOLKIT_PARSE_HPP_INCLUDED
