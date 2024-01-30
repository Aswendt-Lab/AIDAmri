'''
Created on 20.08.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Read Bruker ParaVision JCAMP parameter files (e.g. acqp, method, visu_pars).
'''

from __future__ import print_function

VERSION = 'pv_parser.py v 1.0.2 20200820'

import re
import sys

import collections

import numpy as np

def strfind(string, sub):
    len_sub = len(sub)
    result = []
    if (len_sub == 0) or (len_sub > len(string)):
        return result
    pos = string.find(sub)
    while pos >= 0:
        result.append(pos)
        pos = string.find(sub, pos + len_sub)

    return result

def strtok(string, delimiters=None):
    token = ''
    remainder = ''

    len_str = len(string)
    if len_str == 0:
        return (token, remainder)

    if delimiters is None: # whitespace characters
        delimiters = list(map(chr, list(range(9, 14)) + [32]))

    i = 0
    while string[i] in delimiters:
        i += 1
        if i >= len_str:
            return (token, remainder)

    start = i
    while string[i] not in delimiters:
        i += 1
        if i >= len_str:
            break

    token = string[start:i]
    remainder = string[i:len_str]

    return (token, remainder)

def extract_jcamp_strings(string, get_all=True):
    if string is None:
        result = None
    elif get_all:
        result = re.findall(r'<(.*?)>', string)
    else:
        result = re.search(r'<(.*?)>', string)
        if result is not None:
            result = result.group(1)

    return result

def extract_unit_string(string):
    if string is None:
        result = None
    else:
        result = re.search(r'\[(.*?)\]', string)
        if result is not None:
            result = result.group(1)
        else:
            result = string

    return result

def replace_jcamp_strings(string):
    pos_stop = 0
    elements = []
    str_list = []
    index = 0
    while True:
        pos_start = string.find('<', pos_stop)
        if pos_start < 0:
            elements.append(string[pos_stop:])
            break
        elements.append(string[pos_stop:pos_start])
        pos_stop = string.find('>', pos_start + 1)
        if pos_stop < 0:
            elements.append(string[pos_start:])
            break
        pos_stop += 1
        elements.append(''.join(['<#', str(index), '>']))
        str_list.append(string[pos_start:pos_stop])
        index += 1

    return (''.join(elements), str_list)

def check_struct_list(values, str_list):
    flag_int = True
    flag_float = True
    for value in values:
        if flag_int:
            try:
                value = int(value)
            except ValueError:
                flag_int = False
            else:
                continue
        try:
            value = float(value)
        except ValueError:
            flag_float = False
            break
    if flag_int:
        return (list(map(int, values)), 0)
    if flag_float:
        return (list(map(float, values)), 0)
    # Restore JCAMP strings
    count = len(str_list)
    if count > 0:
        for index, value in enumerate(values):
            result = re.findall(r'<#(.*?)>', value)
            if len(result) == 1:
                str_id = int(result[0])
                values[index] = str_list[str_id]
                count -= 1
                if count == 0:
                    break
            elif len(result) > 1:
                sys.exit("Found more than one ID string in a value: %s" % (value,))

    return (values, len(str_list) - count)

def create_struct_list(string, str_list, restored):
    if len(string) < 1:
        return ([], restored)
    # Split one struct in its parts
    #items = re.split(r'^ +| *, *| +$', string)
    items = re.split(r'(?:^ +| *),(?: *| +$)', string)
    #items = [x.strip(' ') for x in string.split(',')]
    for index, item in enumerate(items):
        #values = re.findall(r'[^\s]+', item)
        values = item.split(' ')
        #values = item.split()
        values, number = check_struct_list(values, str_list)
        if len(values) == 1:
            items[index] = values[0]
        else:
            items[index] = values
        restored += number

    return (items, restored)

def push_list(level, obj_list, obj):
    while level > 0:
        obj_list = obj_list[-1]
        level -= 1
    obj_list.append(obj)

def parse_struct(string, str_list):
    level = 0
    restored = 0
    obj_list = []
    pos_start = string.find('(')
    if pos_start < 0:
        return (obj_list, restored)
    pos_left, start_left = (pos_start + 1, True)
    pos_start = string.find('(', pos_left)
    pos_stop = string.find(')', pos_left)
    while True:
        if (pos_start >= pos_left) and (pos_stop >= pos_left):
            pos_right, start_right = (pos_start, True) if pos_start < pos_stop else (pos_stop, False)
        elif pos_start >= pos_left:
            pos_right, start_right = (pos_start, True)
        elif pos_stop >= pos_left:
            pos_right, start_right = (pos_stop, False)
        else:
            pos_right, start_right = (len(string), False)

        sub = string[pos_left:pos_right].strip(' ')
        if sub.startswith(','):
            sub = sub[1:].lstrip(' ')
        if sub.endswith(','):
            sub = sub[:-1].rstrip(' ')
        #print("sub:%d:%s:" % (len(sub), sub))
        items, restored = create_struct_list(sub, str_list, restored)
        if start_left:
            push_list(level, obj_list, items)
            if start_right:
                level += 1
        else:
            for item in items:
                push_list(level, obj_list, item)
            if not start_right:
                level -= 1
        if pos_right >= len(string):
            break

        pos_left, start_left = (pos_right + 1, start_right)
        if start_left:
            pos_start = string.find('(', pos_left)
        else:
            pos_stop = string.find(')', pos_left)

    return (obj_list, restored)

def check_array_list(values):
    flag_int = True
    flag_float = True
    for value in values:
        if flag_int:
            try:
                value = int(value)
            except ValueError:
                flag_int = False
            else:
                continue
        try:
            value = float(value)
        except ValueError:
            flag_float = False
            break
    if flag_int:
        return np.array(values, dtype=np.int32)
    if flag_float:
        return np.array(values, dtype=np.float64)
    return np.array(values, dtype=object)

def get_array_values(label, sizes, data):
    # Removing whitespaces at the edge of strings
    #data = data.replace('< ', '<')
    #data = data.replace(' >', '>')
    if data.startswith('<'): # Checking if array is a single string or an array of strings ...
        #data = data.replace('> <', '><')
        #values = re.findall(r'<(.*?)>', data)
        values = re.findall(r'<.*?>', data)
        if len(sizes) > 1:
            values = np.array(values, dtype=object)
            if np.prod(sizes[:-1]) == values.size:
                values = values.reshape(sizes[:-1])
        elif len(values) == 1:
            values = values[0]
    elif data.startswith('('): # ... or a struct or an array of structs ...
        if len(sizes) > 1:
            print("Warning: The sizes dimension is greater than 1 for the %s array of structs." % (label,), file=sys.stderr)
        data, str_list = replace_jcamp_strings(data)
        values, restored = parse_struct(data, str_list)
        if len(str_list) != restored:
            print("%s:" % (label,), values)
            sys.exit("Not all replaced JCAMP strings are restored (%d of %d)." % (restored, len(str_list)))
    else: # ... or a simple array (most frequently numeric)
        values = re.findall(r'[^\s]+', data)
        #values = data.split()
        values = np.reshape(check_array_list(values), sizes)

    return values

def read_param_file(filename):
    # Open parameter file
    try:
        fid = open(filename, 'r')
    except IOError as V:
        if V.errno == 2:
            sys.exit("Cannot open parameter file %s" % (filename,))
        else:
            raise

    # Generate header information
    header = collections.OrderedDict()
    weekdays = ('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')

    line = ''
    for index, line in enumerate(fid):
        line = line.lstrip(' \t').rstrip('\r\n')
        if line.startswith('##$'):
            break
        #print("line:%d:%s:" % (len(line), line))

        if line.startswith('##'): # It's a variable with ##
            # Retrieve the Labeled Data Record
            label, value = strtok(line, delimiters='=')
            label = strtok(label, delimiters='#')[0].strip()
            value = strtok(value, delimiters='=')[0].strip()
            # Save value without $
            #value = strtok(value, delimiters='$')[0].strip()
            header[label] = value
        elif line.startswith('$$'): # It's a comment
            comment = strtok(line, delimiters='$')[0].strip()
            if comment.startswith('/'):
                header['Path'] = comment
            elif comment.startswith('process'):
                header['Process'] = comment[8:]
            else:
                pos = strfind(comment[:10], '-')
                if (comment[:3] in weekdays) or ((comment[:2] in ('19', '20')) and (len(pos) == 2)):
                    header['Date'] = comment
                else:
                    header['Header' + str(index + 1)] = comment

    # Check if using a supported version of JCAMP file format
    if 'JCAMPDX' in header:
        version = float(header['JCAMPDX'])
    elif 'JCAMP-DX' in header:
        version = float(header['JCAMP-DX'])
    else:
        sys.exit("The file header is not correct.")

    if (version != 4.24) and (version != 5):
        print("Warning: JCAMP version %s is not supported (%s)." % (version, filename), file=sys.stderr)

    params = collections.OrderedDict()

    # Loop for reading parameters
    while line.lstrip(' \t').startswith('##'):
        result = re.search(r'##(.*)=(.*)', line)
        result = [] if result is None else list(result.groups())

        # Checking if label present and removing proprietary tag
        try:
            label = result[0]
        except:
            label = None
        else:
            if label.startswith('$'):
                label = label[1:]
            #print("label:%d:%s:" % (len(label), label))

        # Checking if value present otherwise value is set to empty string
        try:
            value = result[1]
        except:
            value = ''
        #print("value:%d:%s:" % (len(value), value))

        flag_comment = True if '$$' in line else False

        line = ''
        data = []
        for line in fid:
            if line.lstrip(' \t').startswith('##'):
                break
            if not line.lstrip(' \t').startswith('$$'): # Skip comment line
                if (not flag_comment) and ('$$' in line):
                    flag_comment = True
                #data.append(line.rstrip('\\\r\n'))
                data.append(line.rstrip('\r\n'))
                #print("line:%d:%s:" % (len(data[-1]), data[-1]))

        # Create data string
        data = ''.join(data)
        #print("data:%d:%s:" % (len(data), data))

        if flag_comment:
            sys.exit("Found JCAMP comment ('$$') in LDR %s." % (label,))

        # Checking for END tag
        if (label is None) or (label == 'END'):
            break

        # Checking if value is a string or an array, a struct or a single value
        if value.startswith('( <'):
            print("Warning: The parsing of the LDR %s failed." % (label,), file=sys.stderr)
        elif value.startswith('( '): # A single string, an array of strings or structs or a simple array
            sizes = [int(x) for x in value.strip('( )').split(',')]
            params[label] = get_array_values(label, sizes, data)
        elif value.startswith('('): # A struct
            data = ''.join([value, data])
            params[label] = get_array_values(label, [1], data)[0]
        else: # A single value
            try:
                params[label] = int(value)
            except ValueError:
                try:
                    params[label] = float(value)
                except ValueError:
                    params[label] = value

    fid.close()

    if label != 'END':
        sys.exit("Unexpected end of file: Missing END Statement")

    return (header, params)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Read ParaVision parameter file')
    parser.add_argument('filename', help='ParaVision parameter file (acqp, method, visu_pars)')
    args = parser.parse_args()

    # read parameter file
    header, params = read_param_file(args.filename)

    for (label, value) in header.items():
        print("%s: %s" % (label, value))

    for (label, value) in params.items():
        if isinstance(value, np.ndarray):
            print("%s:" % (label,))
            print(value)
        else:
            print("%s: %s" % (label, value))

if __name__ == '__main__':
    main()
