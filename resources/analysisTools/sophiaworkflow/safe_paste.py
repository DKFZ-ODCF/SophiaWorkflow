#!/usr/bin/env python
#
# Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Similar to coreutils' `paste` but checks that all files indeed have the same number of lines.
#
# Usage:
#
# safe_paste.py file {file}*
#
import sys


class InputError(Exception):

    def __init__(self, message: str):
        self.message = message


def assert_same_linenumber(chunks: list):
    file_number = len(chunks)
    for file_i in range(file_number):
        chunk = chunks[file_i]
        if len(chunk) != len(chunks[0]):
            raise InputError("Inconsistent number of lines: '{}' vs. '{}'".format(files[0], files[file_i]))


def print_chunks(chunks: list, delimiter):
    file_number = len(chunks)
    chunk_size = len(chunks[0])
    for line_offset in range(chunk_size):
        result_line = ""
        for file_i in range(file_number):
            result_line += chunks[file_i][line_offset].rstrip()
            if file_i != file_number - 1:
                result_line += delimiter
        print(result_line)


def all_chunks_empty(chunks: list):
    return all(list(map(lambda chunk: len(chunk) == 0, chunks)))


def any_chunk_not_full(chunks: list, max_chunk_size):
    return any(list(map(lambda chunk: len(chunks) < max_chunk_size, chunks)))


files = sys.argv[1:len(sys.argv)]

try:
    if len(files) == 0:
        raise InputError("Empty argument list. Usage: 'safe_paste.py file {file}*")

    max_chunk_size = 1000
    delimiter = "\t"

    file_handles = list(map(lambda f: open(f, "r"), files))

    while True:
        chunks = list(map(lambda fh: fh.readlines(max_chunk_size), file_handles))

        # For the special case that all files are empty.
        if all_chunks_empty(chunks):
            break

        # Otherwise detect inconsistent line-numbers
        assert_same_linenumber(chunks)
        print_chunks(chunks, delimiter)

        # Stop if the last chunk was read (they have all the same non-null lengths according to previous checks)
        if any_chunk_not_full(chunks, max_chunk_size):
            break

except InputError as e:
    print(e.message, file=sys.stderr)
    exit(1)
