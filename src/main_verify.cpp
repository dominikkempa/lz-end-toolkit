/**
 * @file    src/main_verify.cpp
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
#include <getopt.h>
#include <unistd.h>

#include "lz_end_toolkit_src/verify.hpp"


char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... PARSEFILE TEXTFILE\n"
"Verify if the LZ-End parsing stored in PARSEFILE correctly encodes TEXTFILE\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n",

    program_name);

  std::exit(status);
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",     no_argument,       NULL, 'h'},
    {NULL,       0,                 NULL, 0}
  };

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "h",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(EXIT_FAILURE);
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: PARSEFILE and TEXTFILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the parsing filename.
  std::string parsing_filename = std::string(argv[optind++]);

  if (optind >= argc) {
    fprintf(stderr, "Error: TEXTFILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);

  if (optind < argc) {
    fprintf(stderr, "Error: too many files provided!\n");
    usage(EXIT_FAILURE);
  }

  // Check for the existence of text.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: text file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }

  // Check for the existence of parsing.
  if (!file_exists(parsing_filename)) {
    fprintf(stderr, "Error: parsing file (%s) does not exist\n\n",
        parsing_filename.c_str());
    usage(EXIT_FAILURE);
  }

  verify(parsing_filename, text_filename);
}
