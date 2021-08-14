LZ-End Toolkit
==============


Description
-----------

This package contains implementation of the external-memory
algorithm to construct the LZ-End parsing. The algorithm is
described in the paper

  Dominik Kempa and Dmitry Kosolobov:
  LZ-End Parsing in Compressed Space.
  Data Compression Conference (DCC), IEEE, 2017.

The latest version of the algorithm is available from
https://github.com/dominikkempa/lz-end-toolkit



Compilation and usage
---------------------

The package contains a single Makefile in the main directory.
Type 'make' to build all executables. Run each of the programs
without any argument to see usage instructions.

### Example

Assume the input file is /data/input.dat. The following
demonstrates the simplest usage of the parsing algorithm:

    $ ./parse /data/input.dat
    $ ./verify /data/input.dat.lzend /data/input.dat

The first command computes the LZ-End parsing of the input.
The second command verifies if the computed parsing indeed
correctly encodes the input text (recall, that the parser
is a Monte-Carlo algorithm). If the verifier result is
negative (which can be proved to be very unlikely), the
commands should be repeated. By default, the parsing is
stored in the file named as input with the appended ".lzend"
suffix, the default integer type used to encode phrases
is 40-bit wide (thus enabling the parsing of files of size
up to 1TiB), and the default limit on the phrase length is
2^20. A more advanced usage is demonstrated below.

    $ ./parse /data/input.dat -i 6 -l 4mi -v -o /data/parsing.lzend
    $ ./verify /data/parsing.lzend /data/input.dat

Explanation:

- The -i flag allows specifying the integer size (in bytes)
  used to encode the output parsing. In this example, the
  type is set to 6-byte integer, enabling the parsing of
  inputs up to 256TiB. Currently supported are values from
  the range [4, 8].
- The -l flag allows specifying the limit on the phrase
  length (see the paper above for details on how the limit
  affects the parsing; the default value is sufficient for
  essentially all applications). In this example, the limit
  is set to 4mi = 4 * 2^20.
- the -o flag allows specifying the location and name of
  the file with the output parsing.
- the -v flag enables the verbose mode (i.e., more detailed
  messages during the computation).

Notes:

- The argument of the -l flag (limit on the phrase length)
  can be specified either explicitly or using common suffixes
  such as K, M, G, T, Ki, Mi, Gi, Ti which correspond to
  multipliers: 10^3, 10^6, 10^9, 10^12, 2^10, 2^20 2^30, 2^40.
  Suffix names are not case-sensitive, e.g., Ti = ti, k = K.
- The above flags specifying integer type, output filename,
  etc. can be given in any order.
- Filenames passed as arguments to all programs can be given
  as absolute, relative, and common (such as $HOME) paths,
  e.g., ../input.txt and ~/data/input.txt are valid paths.



Parsing format and decoding
---------------------------

The following describes the format of the parsing. The
format has not been optimized for space, but rather for
the simplicity of decoding.

The parsing file starts with a 8-byte header. The two
least significant bytes of the header store (starting with
the least significant byte): the number of bits used to
encode text symbols minus 1, and the number of bits used
to encode integers in the parsing minus 1. Following the
header is the sequence of triples encoding the parsing.
Each triple consist of a single symbol CHAR and two unsigned
integers (the sizes of all objects are as specified in the
header), in this order. The first of the two integers
stores the ID of the previous phrase and the second integer
stores the length LEN of the phrase. Assume first that the
decoding of the parsing is performed left-to-right and that
all the phrases to the left of the current one has been
decoded. To decode the next phrase (CHAR, ID, LEN) the
decoder first checks if LEN is 1.

- If yes, it appends CHAR at the end of text and proceeds
  with the next phrase.
- If LEN is > 1, The decoder first needs to copy a substring
  of length LEN - 1 ending at the phrase with number ID
  (phrases are numbered left-to-right, starting from 0)
  and append at the end of decoded text. Then, the decoder
  appends CHAR at the end of text and proceeds with the next
  phrase. The end of a phrase is understood as the text
  position containing the CHAR symbol of the phrase.

The simplest implementation of the decoder explicitly
stores the decoded text and the mapping from phrases to
text positions (to determine the end of a phrase a given ID).
This requires the text to be stored in RAM during decoding,
but the parsing can be streamed from disk using a small
buffer. A more space-efficient decoder can copy the LEN - 1
symbols directly from the parsing by using its recursive
structure. This requires the parsing to be kept in RAM, but
allows the output text to be streamed directly to disk.

For cases where LZ-End is used as a compressor as well as
to demonstrate how to read and decode the parsing, this
package provides a simple LZ-End decoder. It implements
a more space-efficient variant of the algorithm mentioned
above. We refer to ./include/lz_end_toolkit/decode.hpp
for details.

### Example

To decode the parsing stored in /data/parsing.lzend, type:

    $ ./decode /data/parsing.lzend

By default, the decoded text is stored in the file named
as the parsing file with the appended ".decoded" suffix.
Other location and name can be specified with the -o flag
(analogous to the -o flag of the parser, see above).



Limitations
-----------

- The current code supports only inputs over byte alphabet.
  The algorithm described in the paper, however, works for
  arbitrary alphabets. Future releases of this package
  will most likely support large alphabet inputs.



Terms of use
------------

If you use this code, please cite the paper mentioned above.
LZ-End Toolkit is released under the MIT/X11 license. See
the file LICENCE for more details.



Authors
-------

LZ-End Toolkit was implemented by:

- [Dominik Kempa](https://scholar.google.com/citations?user=r0Kn9IUAAAAJ)
- [Dmitry Kosolobov](https://scholar.google.com/citations?user=L5boL7MAAAAJ)

