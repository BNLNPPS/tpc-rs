#------------------------------------------------------------------------------
# pycparser: c_json.py
#
# by Michael White (@mypalmike)
#
# This example includes functions to serialize and deserialize an ast
# to and from json format. Serializing involves walking the ast and converting
# each node from a python Node object into a python dict. Deserializing
# involves the opposite conversion, walking the tree formed by the
# dict and converting each dict into the specific Node object it represents.
# The dict itself is serialized and deserialized using the python json module.
#
# The dict representation is a fairly direct transformation of the object
# attributes. Each node in the dict gets one metadata field referring to the
# specific node class name, _nodetype. Each local attribute (i.e. not linking
# to child nodes) has a string value or array of string values. Each child
# attribute is either another dict or an array of dicts, exactly as in the
# Node object representation. The "coord" attribute, representing the
# node's location within the source code, is serialized/deserialized from
# a Coord object into a string of the format "filename:line[:column]".
#
# Example TypeDecl node, with IdentifierType child node, represented as a dict:
#     "type": {
#         "_nodetype": "TypeDecl",
#         "coord": "c_files/funky.c:8",
#         "declname": "o",
#         "quals": [],
#         "type": {
#             "_nodetype": "IdentifierType",
#             "coord": "c_files/funky.c:8",
#             "names": [
#                 "char"
#             ]
#         }
#     }
#------------------------------------------------------------------------------
from __future__ import print_function

import json
import sys
import re

# This is not required if you've installed pycparser into
# your site-packages/ with setup.py
#
sys.path.extend(['.', '..'])

from pycparser import parse_file, c_ast
from pycparser.plyparser import Coord


RE_CHILD_ARRAY = re.compile(r'(.*)\[(.*)\]')
RE_INTERNAL_ATTR = re.compile('__.*__')


class CJsonError(Exception):
    pass


def memodict(fn):
    """ Fast memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = fn(key)
            return ret
    return memodict().__getitem__


@memodict
def child_attrs_of(klass):
    """
    Given a Node class, get a set of child attrs.
    Memoized to avoid highly repetitive string manipulation

    """
    non_child_attrs = set(klass.attr_names)
    all_attrs = set([i for i in klass.__slots__ if not RE_INTERNAL_ATTR.match(i)])
    return all_attrs - non_child_attrs


def to_dict(node):
    """ Recursively convert an ast into dict representation. """
    klass = node.__class__

    result = {}

    # Metadata
    result['_nodetype'] = klass.__name__

    # Local node attributes
    for attr in klass.attr_names:
        result[attr] = getattr(node, attr)

    # Coord object
    if node.coord:
        result['coord'] = str(node.coord)
    else:
        result['coord'] = None

    # Child attributes
    for child_name, child in node.children():
        # Child strings are either simple (e.g. 'value') or arrays (e.g. 'block_items[1]')
        match = RE_CHILD_ARRAY.match(child_name)
        if match:
            array_name, array_index = match.groups()
            array_index = int(array_index)
            # arrays come in order, so we verify and append.
            result[array_name] = result.get(array_name, [])
            if array_index != len(result[array_name]):
                raise CJsonError('Internal ast error. Array {} out of order. '
                    'Expected index {}, got {}'.format(
                    array_name, len(result[array_name]), array_index))
            result[array_name].append(to_dict(child))
        else:
            result[child_name] = to_dict(child)

    # Any child attributes that were missing need "None" values in the json.
    for child_attr in child_attrs_of(klass):
        if child_attr not in result:
            result[child_attr] = None

    return result


def to_json(node, **kwargs):
    """ Convert ast node to json string """
    return json.dumps(to_dict(node), **kwargs)


def file_to_dict(filename):
    """ Load C file into dict representation of ast """
    ast = parse_file(filename, use_cpp=True)
    return to_dict(ast)


def file_to_json(filename, **kwargs):
    """ Load C file into json string representation of ast """
    ast = parse_file(filename, use_cpp=True)
    return to_json(ast, **kwargs)


def _parse_coord(coord_str):
    """ Parse coord string (file:line[:column]) into Coord object. """
    if coord_str is None:
        return None

    vals = coord_str.split(':')
    vals.extend([None] * 3)
    filename, line, column = vals[:3]
    return Coord(filename, line, column)


def _convert_to_obj(value):
    """
    Convert an object in the dict representation into an object.
    Note: Mutually recursive with from_dict.

    """
    value_type = type(value)
    if value_type == dict:
        return from_dict(value)
    elif value_type == list:
        return [_convert_to_obj(item) for item in value]
    else:
        # String
        return value


def from_dict(node_dict):
    """ Recursively build an ast from dict representation """
    class_name = node_dict.pop('_nodetype')

    klass = getattr(c_ast, class_name)

    # Create a new dict containing the key-value pairs which we can pass
    # to node constructors.
    objs = {}
    for key, value in node_dict.items():
        if key == 'coord':
            objs[key] = _parse_coord(value)
        else:
            objs[key] = _convert_to_obj(value)

    # Use keyword parameters, which works thanks to beautifully consistent
    # ast Node initializers.
    return klass(**objs)


def from_json(ast_json):
    """ Build an ast from json string representation """
    return from_dict(json.loads(ast_json))


def extract_decls(ast_dict, decl=None):
    """Extracts struct member declarations defined in a header file. A single
    struct definition is assumed per file
    """

    class Decl:
        """A trivial class with primary purpose to provide access to data attributes
        by dot notation used in the `extract_decls` function
        """
        def __init__(self):
            self.name = None
            self.type = None
            self.dims = []

    if isinstance(ast_dict, dict):
        for key, value in ast_dict.items():
            if key is '_nodetype' and value is "Decl":
                decl = Decl()
            elif key is 'init':
                decl.dims = decl.dims[::-1]
                yield decl
            elif key is 'declname':
                decl.name = value
            elif key is 'dim':
                decl.dims += [int(value['value'])]
            elif key is 'names':
                decl.type = ' '.join(value)
            elif isinstance(value, dict):
                for d in extract_decls(value, decl):
                    yield d
            elif isinstance(value, list) or isinstance(value, tuple):
                for v in value:
                    for d in extract_decls(v, decl):
                        yield d


def generate_yamlcpp(ast_dict, pre=None):
    """Generates c++ code for converting C structures to YAML format and back by
    the `yaml-cpp` library
    """

    # Skip header files without the proper struct
    if not ast_dict['ext'] or not len(ast_dict['ext']):
       return ""

    cxx_code_encode = ""
    cxx_code_decode = ""

    for decl in extract_decls(ast_dict):
        if not decl.name or not decl.type: continue

        if decl.dims:
            from functools import reduce
            size = reduce(lambda x, y: x*y, decl.dims)
            if decl.type == "unsigned char" and len(decl.dims) == 1:
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"] = vector<int>(st.{decl.name}, st.{decl.name} + {size});"
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"].SetStyle(YAML::EmitterStyle::Flow);"
                cxx_code_decode += f"\n\tauto {decl.name} = node[\"{decl.name}\"].as< vector<int> >();"
                cxx_code_decode += f"\n\tcopy(begin({decl.name}), end({decl.name}), st.{decl.name});"
            elif decl.type == "char" and len(decl.dims) == 1:
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"] = string(st.{decl.name});"
                cxx_code_decode += f"\n\tauto {decl.name} = node[\"{decl.name}\"].as<string>();"
                cxx_code_decode += f"\n\tcopy(begin({decl.name}), end({decl.name}), st.{decl.name});"
                cxx_code_decode += f"\n\tst.{decl.name}[{decl.name}.size()] ='\\0';"
            else:
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"] = reinterpret_cast<const array<{decl.type}, {size}>&>( st.{decl.name} );"
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"].SetStyle(YAML::EmitterStyle::Flow);"
                cxx_code_decode += f"\n\tauto {decl.name} = node[\"{decl.name}\"].as<array<{decl.type}, {size}>>();"
                cxx_code_decode += f"\n\tcopy(begin({decl.name}), end({decl.name}), reinterpret_cast<{decl.type}*>(st.{decl.name}));"
        else:
            if decl.type == "unsigned char":
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"] = static_cast<int>(st.{decl.name});"
                cxx_code_decode += f"\n\tst.{decl.name} = static_cast<{decl.type}>(node[\"{decl.name}\"].as<int>());"
            else:
                cxx_code_encode += f"\n\tnode[\"{decl.name}\"] = st.{decl.name};"
                cxx_code_decode += f"\n\tst.{decl.name} = node[\"{decl.name}\"].as<{decl.type}>();"

    struct_name = ast_dict['ext'][0]['type']['name']
    # Pop the last three charasters "_st"
    header_name = struct_name[:-3]

    cxx_code_common = """
    #include <algorithm>
    #include <cstdint>
    #include "yaml-cpp/yaml.h"

    namespace YAML {
    using std::copy;
    using std::begin;
    using std::end;
    using std::array;
    using std::string;
    using std::vector;
    }
    """

    cxx_code = """
    #include "%(header_name)s.h"

    namespace YAML {
    template<>
    struct convert<%(struct_name)s> {
      static Node encode(const %(struct_name)s& st) {
        Node node;
        %(cxx_code_encode)s
        return node;
      };

      static bool decode(const Node& node, %(struct_name)s& st) {
        if(!node.IsMap()) {
          return false;
        }
        %(cxx_code_decode)s
        return true;
      }
    };
    }
    """ % {'header_name': header_name, 'struct_name': struct_name,
           'cxx_code_encode': cxx_code_encode, 'cxx_code_decode': cxx_code_decode}

    import inspect
    return inspect.cleandoc(cxx_code)


#------------------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Some test code...
        # Do trip from C -> ast -> dict -> ast -> json, then print.
        ast_dict = file_to_dict(sys.argv[1])
        cxx_code = generate_yamlcpp(ast_dict)
        print(cxx_code)
    else:
        print("Please provide a filename as argument")
