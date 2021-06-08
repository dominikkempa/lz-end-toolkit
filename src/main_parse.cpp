/**
 * @file    src/main_parse.cpp
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <cctype>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>

#include "lz_end_toolkit_src/parse.hpp"
#include "uint24.hpp"
#include "uint40.hpp"
#include "uint48.hpp"
#include "uint56.hpp"


char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the LZ-End parsing of text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes. Int type needs to be wide\n"
"                          enough to encode any position in FILE. Default: 5.\n"
"                          Currently supported values: 4-8\n"
"  -l, --limit=LEN         set limit for phrase length to LEN. Metric and IEC\n"
"                          suffixes are recognized, e.g., -l 10k, -l 1Mi, -l 3G\n"
"                          gives LEN = 10^4, 2^20, 3*10^6. Default: 1Mi\n"
"  -o, --output=OUTFILE    specify output filename. Default: FILE.lzend\n"
"  -v, --verbose           enable verbose progress messages during computation\n",

    program_name);

  std::exit(status);
}

bool file_exists(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

template<typename int_type>
bool parse_number(char *str, int_type *ret) {
  *ret = 0;
  std::uint64_t n_digits = 0;
  std::uint64_t str_len = std::strlen(str);
  while (n_digits < str_len && std::isdigit(str[n_digits])) {
    std::uint64_t digit = str[n_digits] - '0';
    *ret = (*ret) * 10 + digit;
    ++n_digits;
  }

  if (n_digits == 0)
    return false;

  std::uint64_t suffix_length = str_len - n_digits;
  if (suffix_length > 0) {
    if (suffix_length > 2)
      return false;

    for (std::uint64_t j = 0; j < suffix_length; ++j)
      str[n_digits + j] = std::tolower(str[n_digits + j]);
    if (suffix_length == 2 && str[n_digits + 1] != 'i')
      return false;

    switch(str[n_digits]) {
      case 'k':
        if (suffix_length == 1)
          *ret *= 1000;
        else
          *ret <<= 10;
        break;
      case 'm':
        if (suffix_length == 1)
          *ret *= 1000000;
        else
          *ret <<= 20;
        break;
      case 'g':
        if (suffix_length == 1)
          *ret *= 1000000000;
        else
          *ret <<= 30;
        break;
      case 't':
        if (suffix_length == 1)
          *ret *= 1000000000000;
        else
          *ret <<= 40;
        break;
      default:
        return false;
    }
  }

  return true;
}

template<typename char_type>
void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose,
    std::uint64_t int_size) {
  if (int_size == 4)
    parse<char_type, std::uint32_t>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else if (int_size == 5)
    parse<char_type, uint40>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else if (int_size == 6)
    parse<char_type, uint48>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else if (int_size == 7)
    parse<char_type, uint56>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else if (int_size == 8)
    parse<char_type, std::uint64_t>(text_filename,
        output_filename, phrase_length_limit, verbose);
  else {
    fprintf(stderr, "\nError: parse: unsupported int type!\n");
    std::exit(EXIT_FAILURE);
  }
}

void parse(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t phrase_length_limit,
    bool verbose,
    std::uint64_t char_size,
    std::uint64_t int_size) {
  if (char_size == 1)
    parse<std::uint8_t>(text_filename, output_filename,
        phrase_length_limit, verbose, int_size);
  else {
    fprintf(stderr, "\nError: parse: unsupported char type!\n");
    std::exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",     no_argument,       NULL, 'h'},
    {"intsize",  required_argument, NULL, 'i'},
    {"limit",    required_argument, NULL, 'l'},
    {"output",   required_argument, NULL, 'o'},
    {"verbose",  no_argument,       NULL, 'v'},
    {NULL,       0,                 NULL, 0}
  };

  std::uint64_t char_size = 1;
  std::uint64_t int_size = 5;
  std::uint64_t phrase_length_limit = ((std::uint64_t)1 << 20);
  std::string output_filename("");
  bool verbose = false;

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "hi:l:o:v",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(EXIT_FAILURE);
      case 'i':
        int_size = std::atol(optarg);
        if (int_size < 4 || int_size > 8) {
          fprintf(stderr, "Error: invalid int size (%lu)\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'l':
        {
          bool ok = parse_number(optarg, &phrase_length_limit);
          if (!ok) {
            fprintf(stderr, "Error: parsing phrase length "
                "limit (%s) failed\n\n", optarg);
            usage(EXIT_FAILURE);
          }
          break;
        }
      case 'o':
        output_filename = std::string(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  if (output_filename.empty())
    output_filename = text_filename + ".lzend";

  // Check for the existence of text.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }

  if (file_exists(output_filename)) {

    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          output_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }

  // Run the algorithm.
  parse(text_filename, output_filename,
      phrase_length_limit, verbose, char_size, int_size);
}
