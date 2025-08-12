# A collection of utility methods and classes that can be useful in other projects as well
from __future__ import annotations
import io
import re
import shutil
import zipfile
from collections import namedtuple
from collections.abc import Callable
from pathlib import Path
import argparse
import os

class ValidationError(Exception):
    pass

class DiscontinuousRange:
    def __init__(self, *items: int | range, range_str: str = None, min_value: int = None, max_value: int = None):
        """
            Creates a discontinous range generator. Items are given in the order that they are inputted into the range. Use either *items or range_str
            :param items: a collection of elements to turn into a range
            :param range_str: a string composed of elements separated by commas. The elements can be integers or ranges (denoted with a -, end is INCLUSIVE)
            :param min_value: for user-input, the minimum allowed value in the DiscontinuousRange
            :param max_value: for user-input, the maximum allowed value in the DiscontinuousRange (exclusive)
            :return:
        """
        if len(items) == 0: items = range_str.replace(" ", "").split(',')
        self.items = [self._get_and_verify(item, min_value=min_value, max_value=max_value) for item in items]

    @staticmethod
    def _get_and_verify(item: str | int | range, min_value: int = None, max_value: int = None):
        range_obj = _get_range(item)
        first_elem, last_elem = (range_obj.start, range_obj.stop - 1) if range_obj.step > 0 else (range_obj.stop + 1, range_obj.start)
        if is_not_below_min(min_value, first_elem) and is_below_max(max_value, last_elem):
            return range_obj
        raise ValueError(f"Range {range_obj} is not between {min_value or 'negative infinity'} inclusive and {max_value or 'infinity'} exclusive")

    @staticmethod
    def template(min_value = None, max_value = None, input_string = True):
        if input_string: return lambda x: DiscontinuousRange(range_str=x, min_value=min_value, max_value=max_value)
        return lambda *x: DiscontinuousRange(*x, min_value=min_value, max_value=max_value)

    def __iter__(self):
        return (i for item in self.items for i in item)

    def __reversed__(self):
        return (i for item in reversed(self.items) for i in reversed(item))
    def __str__(self):
        return ", ".join([(str(range.start) if range.start == range.stop-1 else f"{range.start}-{range.stop - 1}") for range in self.items])

    def __repr__(self):
        return f"DiscontinuousRange({str(self)})"

def identity(x):
    return x

def input_value(msg : str, mapper: type[any] = identity, predicate = lambda n: True, fail_message : str = None, initial_value = None, retry_if_fail=True, map_initial_value = False):
    fail_message = fail_message or msg
    while True:
        try:
            if map_initial_value and initial_value is not None: initial_value = mapper(initial_value)
            x = initial_value if initial_value is not None else mapper(input(msg))
            initial_value = None
            if not predicate(x):
                raise ValueError()
            return x
        except Exception as e:
            if retry_if_fail: print(fail_message)
            else: raise ValueError(fail_message) from e

def input_int_in_range(min: int = None, max: int = None, msg : str = None, fail_message : str = None,
                       initial_value = None, extra_predicate=lambda x: False, retry_if_fail=True) -> int:
    fail_message = fail_message or f"Please input an integer between {min} (inclusive) and {max} (exclusive)."
    msg = msg or f"Input an integer between {min} (inclusive) and {max} (exclusive): "
    return input_value(msg, mapper=int, initial_value=initial_value, predicate = lambda x: ((min is None or x >= min) and (max is None or x < max)) or extra_predicate(x),
                       fail_message = fail_message, retry_if_fail=retry_if_fail)

def input_range(min: int = None, max: int = None, msg : str = None, fail_message : str = None,
                       initial_value = None, retry_if_fail=True) -> DiscontinuousRange:
    fail_message = fail_message or f"Please input a range contained between {min} (inclusive) and {max} (exclusive)."
    msg = msg or f"Input a range contained between {min} (inclusive) and {max} (exclusive): "
    return input_value(msg, mapper=DiscontinuousRange.template(min_value=min, max_value=max), initial_value=initial_value,
                       fail_message = fail_message, retry_if_fail=retry_if_fail)

def input_bool(true_vals: str | list[str] = "y,yes,true,t", false_vals: str | list[str] = "n,no,false,f", msg : str = None, fail_message : str = None,
                       initial_value = None, retry_if_fail=True, map_initial_value=True) -> bool:
    fail_message = fail_message or f"Please input a valid boolean (y/n, True/False, etc.)"
    msg = msg or f"Input a valid boolean (y/n, True/False, etc.)"
    true_vals = true_vals if isinstance(true_vals, list) else re.split(",\\s*", true_vals)
    false_vals = false_vals if isinstance(false_vals, list) else re.split(",\\s*", false_vals)
    def matches(string: str):
        string = string.strip().lower()
        if string in true_vals: return True
        if string in false_vals: return False
        raise ValueError(f"String must be one of {true_vals} or {false_vals}")
    return input_value(msg, mapper=matches, initial_value=initial_value,
                       fail_message = fail_message, retry_if_fail=retry_if_fail, predicate=lambda x: isinstance(x, bool), map_initial_value=map_initial_value)

def input_path(*suffixes: str, msg : str = None, fail_message : str = None, initial_value = None, retry_if_fail=True, map_initial_value=False):
    fail_message = fail_message or f"Please input the path to a file. Extension can be: {', '.join(suffixes)}"
    msg = msg or f"Input the path to a file. Extension can be: {', '.join(suffixes)}"
    return input_value(msg, fail_message=fail_message, mapper=Path, predicate=lambda path: path.exists() and path.suffix in suffixes, initial_value=initial_value, retry_if_fail=retry_if_fail, map_initial_value=map_initial_value)

def input_path_string(*suffixes: str, msg : str = None, fail_message : str = None, initial_value = None, retry_if_fail=True, map_initial_value=False):
    fail_message = fail_message or f"Please input the path to a file. Extension can be: {', '.join(suffixes)}"
    msg = msg or f"Input the path to a file. Extension can be: {', '.join(suffixes)}"
    return input_value(msg, fail_message=fail_message, predicate=lambda string: Path(string).exists() and Path(string).suffix in suffixes, initial_value=initial_value, retry_if_fail=retry_if_fail, map_initial_value=map_initial_value)

def remove_files(*files):
    for file in files:
        remove_if_exists(file)

def remove_if_exists(path: Path):
    if path.exists():
        os.remove(path)

def is_empty(directory: Path):
    if not directory.is_dir(): raise ValueError(f"Input must be an existing directory. {str(directory)} either is not a directory or does not exist")
    return not any(directory.iterdir())

def validate_doesnt_throw(func: Callable, *input, msg="", **kwargs):
    try:
        return func(*input, **kwargs)
    except Exception as e:
        raise ValidationError(msg) from e

def validate_arg(boolean: bool, msg):
    if not boolean: raise ValidationError(msg)

def validate_range_arg(input: int, min = None, max = None, name = "input", overwrite_msg = None, extra_predicate = lambda x: False):
    """
    Validate an argument that conforms to a range
    :param input:
    :param min: the max value of the argument, exclusive
    :param max:
    :param msg:
    :return:
    """
    validate_arg(type(input) is int, overwrite_msg or f"The given {name} is not an integer! It must be an integer between {min or 'negative Infinity'} inclusive and {max or 'Infinity'} exclusive")
    min_OK = is_not_below_min(min, input)
    max_OK = is_below_max(max, input)
    validate_arg((max_OK and min_OK) or extra_predicate(input), overwrite_msg or f"The given {name} ({input}) is {'too small' if max_OK else 'too large'}! It must be an integer between {min or 'negative Infinity'} inclusive and {max or 'Infinity'} exclusive")

def is_not_below_min(min, *values):
    return all(min is None or value >= min for value in values)

def is_below_max(max, *values):
    return all(max is None or value < max for value in values)

def optional_argument(request, arg_name: str, cmd_line_name: str = None, default_value = None, type: Callable =None):
    if cmd_line_name is not None:
        if request.form[arg_name] != "":
            return f" {cmd_line_name} {request.form[arg_name]}"
        else:
            return f" {cmd_line_name} {default_value}" if default_value is not None else ""
    elif default_value is not None:
        return default_value if request.form[arg_name] == "" else type(request.form[arg_name])
    else:
        raise ValueError("Must include either cmd_line_name or default_value (or both)")

def get_or_none(obj, attr):
    return getattr(obj, attr, None)

def bounded_int(string, min=None, max=None):
       value = int(string)
       if is_not_below_min(min, value) and is_below_max(max + 1, value):
           return value
       else:
           raise argparse.ArgumentTypeError(
               f'Value not in range {min or "-∞"}-{max or "∞"}. Please either keep it in range or leave it out.')

def path_string(string, suffix=".ct"):
    path = Path(string).resolve()
    if path.exists() and path.suffix in suffix:  # in returns true if string ==
        return str(path)
    else:
        raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

def path_arg(string, suffix=".ct"):
    path = Path(string).resolve()
    if path.exists() and path.suffix in suffix:  # in returns true if string ==
        return path
    else:
        raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

def get_folder_as_zip(folder_path: Path) -> bytes:
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, start=folder_path)
                zipf.write(abs_path, arcname=rel_path)
    zip_buffer.seek(0)
    return zip_buffer.getvalue()

ParsedFile = namedtuple("ParsedFile", ["parent", "stem", "suffix"])
def parse_file_input(filein: str | Path, output_dir = None) -> ParsedFile:
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = output_dir or filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return ParsedFile(mb_userpath, fname, filein_path.suffix)

def safe_remove_tree(folder: Path, root_dir): #just to be really safe when removing files
    if not folder.exists(): return
    if root_dir in folder.parents:
        shutil.rmtree(folder)
    else:
        raise ValueError(f"The given folder {folder} is not a child of root folder {root_dir}")

def _get_range(item: str | range | int) -> range:
    if isinstance(item, range): return item
    if isinstance(item, int): return range(item, item+1)

    pair = item.split('-')  # may only be one item, but that's OK
    return range(int(pair[0]), int(pair[-1])+1)

_Colors = dict(
    RED='\033[91m',
    YELLOW = '\033[93m',
    GREEN='\033[92m',
    BLUE = '\033[94m',
    CYAN = '\033[96m',
    MAGENTA = '\033[95m',
    WHITE = '\033[97m',
    HEADER='\033[95m',
    FAIL = '\033[91m',
    WARNING='\033[93m',
    ENDC = '\033[0m',
    BOLD = '\033[1m',
    UNDERLINE = '\033[4m')

def print_style(msg, *colors):
    print("".join(_Colors[color.upper()] for color in colors) + msg + _Colors['ENDC'])

if __name__ == "__main__":
    print("Debug") #debug