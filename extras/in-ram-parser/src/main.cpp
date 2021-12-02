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
#include <ctime>
#include <cctype>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>

#include "../include/compute_lzend.hpp"
#include "../include/types/uint24.hpp"
#include "../include/types/uint40.hpp"
#include "../include/types/uint48.hpp"
#include "../include/types/uint56.hpp"


char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the LZ-End parsing of text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes (default: 5). Currently\n"
"                          supported values: 4-8\n"
"  -o, --output=OUTFILE    specify output filename. Default: FILE.lzend\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
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

void compute_lzend_for_file(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t int_size) {
  if (int_size == 4)
    compute_lzend_for_file<std::uint8_t, std::uint32_t>(text_filename,
        output_filename);
  else if (int_size == 5)
    compute_lzend_for_file<std::uint8_t, uint40>(text_filename,
        output_filename);
  else if (int_size == 6)
    compute_lzend_for_file<std::uint8_t, uint48>(text_filename,
        output_filename);
  else if (int_size == 7)
    compute_lzend_for_file<std::uint8_t, uint56>(text_filename,
        output_filename);
  else if (int_size == 8)
    compute_lzend_for_file<std::uint8_t, std::uint64_t>(text_filename,
        output_filename);
  else {
    fprintf(stderr, "\nError: compute_lzend_for_file: unsupported int type!\n");
    std::exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",     no_argument,       NULL, 'h'},
    {"intsize",  required_argument, NULL, 'i'},
    {"output",   required_argument, NULL, 'o'},
    {NULL,       0,                 NULL, 0}
  };

  std::uint64_t int_size = 5;
  std::string output_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "hi:o:",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(EXIT_FAILURE);
        break;
      case 'i':
        int_size = std::atol(optarg);
        if (!(int_size == 4 ||
              int_size == 5 ||
              int_size == 6 ||
              int_size == 7 ||
              int_size == 8)) {
          fprintf(stderr, "Error: invalid int size (%lu)\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'o':
        output_filename = std::string(optarg);
        break;
      default:
        usage(EXIT_FAILURE);
        break;
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

  // Run the algorithm.
  compute_lzend_for_file(text_filename, output_filename, int_size);
}

