#!/usr/bin/env python

import sys


def assert_same_linenumber(chunks):
    file_number = len(chunks)
    for file_i in range(file_number):
        chunk = chunks[file_i]
        if len(chunk) != len(chunks[0]):
            raise AssertionError("Inconsistent number of lines: '{}' vs. '{}'".format(files[0], files[file_i]))


def print_chunks(chunks, delimiter):
    file_number = len(chunks)
    chunk_size = len(chunks[0])
    for line_offset in range(chunk_size):
        result_line = ""
        for file_i in range(file_number):
            result_line += chunks[file_i][line_offset].rstrip()
            if file_i != file_number - 1:
                result_line += delimiter
        print(result_line)


def all_chunks_empty(chunks):
    return all(list(map(lambda chunk: len(chunk) == 0, chunks)))


def any_chunk_not_full(chunks, max_chunk_size):
    return any(list(map(lambda chunk: len(chunks) < max_chunk_size, chunks)))


files = sys.argv[1:len(sys.argv)]
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

    # Finally, stop if the last chunk was read (they have all the same non-null lengths according to previous checks)
    if any_chunk_not_full(chunks, max_chunk_size):
        break
