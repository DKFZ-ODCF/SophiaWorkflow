#!/usr/bin/env python
#
# Copyright (C) 2020 Philip R. Kensche, DKFZ Heidelberg
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


def assert_same_linenumber(data: list):
    file_number = len(data)
    for file_i in range(file_number):
        if len(data[file_i]) != len(data[0]):
            raise InputError("Inconsistent number of lines: {} ({}) != {} ({})".
                             format(len(files[0]), files[0],
                                    len(files[file_i]), files[file_i]))


def print_data(all_data: list, delimiter):
    file_number = len(all_data)
    for line_i in range(len(all_data[0])):
        result_line = ""
        for file_i in range(file_number):
            result_line += all_data[file_i][line_i].rstrip()
            if file_i != file_number - 1:
                result_line += delimiter
        print(result_line)


def all_files_empty(all_data: list):
    return all(list(map(lambda file_data: len(file_data) == 0, all_data)))


files = sys.argv[1:len(sys.argv)]

delimiter = "\t"

try:
    if len(files) == 0:
        raise InputError("Empty argument list. Usage: 'safe_paste.py file {file}*")

    file_handles = map(lambda f: open(f, "r"), files)
    all_data = list(map(lambda fh: fh.readlines(), file_handles))

    if not all_files_empty(all_data):
        assert_same_linenumber(all_data)
        print_data(all_data, delimiter)

except InputError as e:
    print(e.message, file=sys.stderr)
    exit(1)
