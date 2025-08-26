# A collection of utility methods and classes specifically for this project
from __future__ import annotations

import shutil
import sys
import zipfile
from argparse import Namespace
from collections import namedtuple
from collections.abc import Callable
from io import UnsupportedOperation
from pathlib import Path
import io
from typing import IO

from . import util
from .util import remove_if_exists, ValidationError, safe_remove_tree, is_empty


def run_command_line(run: Callable, *args, **kwargs):
    try:
        return run(*args, **kwargs)
    except ValidationError as e:
        print(f"{str(e)}. Program terminated.", file=sys.stderr)

class UnclosableStringIO(io.StringIO):
    def close(self):
        # Override close so it doesn't actually close the stream
        pass

class UnclosableBytesIO(io.BytesIO):
    def close(self):
        # Override close so it doesn't actually close the stream
        pass

class FileManager:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.files = dict()
        self.files_to_delete = set()

    def add_file(self, rel_path: str, path: Path, delete_when_clean = True, is_directory = False):
        if delete_when_clean: self.files_to_delete.add(path)
        if not is_directory: self.files[rel_path] = path

    def add_to_delete(self, path: Path):
        self.files_to_delete.add(path)

    def get_files(self) -> dict:
        return self.files

    def remove_file(self, rel_path: str, path: Path):
        remove_if_exists(path)
        self.files_to_delete.discard(path)
        self.files.pop(rel_path, None)

    def cleanup(self):
        directories = set()
        for file in self.files_to_delete:
            if file.is_dir():
                directories.add(file)
            else:
                if file.parent != self.output_dir: directories.add(file.parent)
                remove_if_exists(file)

        for directory in directories:
            if directory.exists() and is_empty(directory): #directory is empty
                safe_remove_tree(directory, self.output_dir)

class ProgramObject:
    def __init__(self, output_dir: Path, file_stem: str, arguments: Namespace, **kwargs):
        self.result_obj = Namespace(**kwargs)
        self.output_dir = output_dir
        self.file_stem = file_stem
        self.arguments = arguments
        self.file_manager = FileManager(output_dir)
        if output_dir is not None: output_dir.mkdir(parents=True, exist_ok=True)

    def save_buffer(self, rel_path: str, register_to_delete=True):
        """
        Returns a buffer to use to store files
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :return:
        """
        self.register_file(rel_path, register_to_delete=register_to_delete)
        return self.file_path(rel_path) #test

    def open_buffer(self, rel_path: str, mode = "w", register_to_delete=True) -> IO:
        """
        Buffer used to mock the open() method
        :param rel_path:
        :param mode:
        :param register_to_delete:
        :return:
        """
        return open(self.save_buffer(rel_path, register_to_delete=register_to_delete), mode)

    def reset_buffer(self, rel_path: str):
        path = self.file_path(rel_path)
        self.file_manager.remove_file(self.format_relative_path(rel_path), path)


    def file_path(self, rel_path: str, register=False, register_to_delete=False, is_directory=False) -> Path:
        """
        Returns a path to use to store files. Same as save_buffer if from the command line, but will always
        return a path (unlike save_buffer which may return a buffer)
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :param register: register this file as a ProgramObject file (to be removed if the program fails)
        :return:
        """
        path = self.output_dir / self.format_relative_path(rel_path)
        if register: self.register_file(rel_path, true_path=path, is_directory=is_directory, register_to_delete=register_to_delete)
        return path

    def format_relative_path(self, rel_path: str) -> str:
        return rel_path.replace("[fname]", self.file_stem)

    def create_dir(self, rel_path: str, register=True):
        path = self.file_path(rel_path, register=register, is_directory=True)
        path.mkdir(parents=True, exist_ok=True)

    def register_file(self, rel_path: str, true_path: Path = None, register_to_delete: bool = False, is_directory: bool=False):
        self.file_manager.add_file(rel_path=self.format_relative_path(rel_path), path=true_path or self.file_path(rel_path),
                                   delete_when_clean=register_to_delete, is_directory=is_directory)
    def register_to_delete(self, rel_path: str = None, true_path: Path = None):
        self.file_manager.add_to_delete(path=true_path or self.file_path(rel_path))

    def get_arg(self, argument):
        return getattr(self.arguments, argument)

    def set_args(self, **kwargs):
        for argument, value in kwargs.items():
            setattr(self.arguments, argument, value)

    def get_result_arg(self, argument):
        return getattr(self.result_obj, argument)

    def set_result_args(self, **kwargs):
        for argument, value in kwargs.items():
            setattr(self.result_obj, argument, value)
    def to_zip(self, name) -> tuple[bytes, str]:
        archive_name = name.replace("[fname]", self.file_stem)
        return util.get_folder_as_zip(self.output_dir), archive_name

    def validate(self, boolean: bool, msg: str):
        if not boolean:
            self.quit_program(msg, error=ValidationError)

    def quit_program(self, msg, error=ValidationError):
        """
        Quit the program. The user should make sure to catch ValidationError when needed in order to prevent an ugly traceback
        (any other error may be treated as an actual error). See the global function run_command_line
        :param msg: the message to print when quitting the program
        :param error: The error to raise in order to quit. Default is ValidationError.
        :return:
        """
        self.cleanup()
        raise error(msg)

    def cleanup(self):
        self.file_manager.cleanup()

#todo: prevent collision in file_dict and buffer_dict
class BufferedProgramObject(ProgramObject):
    def __init__(self, output_dir: Path, file_stem: str, arguments: Namespace, **kwargs):
        ProgramObject.__init__(self, output_dir, file_stem, arguments, **kwargs)
        self.buffer_dict = dict()
    def save_buffer(self, rel_path: str, is_string: bool = True):
        """
        Returns a string buffer to use to store files
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :param is_string: True if using StringIO, false for BytesIO
        :return:
        """
        path = Path(self.format_relative_path(rel_path))
        return self.buffer_dict.setdefault(path, UnclosableStringIO() if is_string else UnclosableBytesIO()) #test

    def open_buffer(self, rel_path: str, mode = "w", register_to_delete=True):
        """

        :param rel_path:
        :param mode: the mode to write in, w to remove the previous content from the buffer. For compatibility with parent class
        :return:
        """
        if "w" in mode: self.reset_buffer(rel_path)
        buffer = self.save_buffer(rel_path, not ('b' in mode))
        if not "a" in mode: buffer.seek(0)
        return buffer

    def reset_buffer(self, rel_path: str):
        path = self.format_relative_path(rel_path)
        if path in self.buffer_dict:
            self.buffer_dict[path].reset()
        if path in self.file_manager.get_files():
            super().reset_buffer(rel_path)

    # def register_file(self, rel_path: str, true_path:  Path = None, is_directory=False, register_to_delete=True):

    def file_path(self, rel_path: str, register=False, register_to_delete=False, is_directory=False):
        """
        Return a file
        :param rel_path:
        :param register:
        :param register_to_delete: if register is true, also register this file to be deleted if the program crashes
        :param type:
        :return:
        """
        if self.output_dir is None:
            raise UnsupportedOperation("Can't get a file path if output_dir is None")
        return super().file_path(rel_path, register=register, is_directory=is_directory, register_to_delete=register_to_delete)

    def to_zip(self, name) -> tuple[bytes, str]:
        archive_name = name.replace("[fname]", self.file_stem)
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for path, content in self.buffer_dict.items():
                zipf.writestr(zinfo_or_arcname=str(path), data=content.getvalue())
            for path, abs_path in self.file_manager.get_files().items():
                zipf.write(filename=abs_path, arcname=path.path)
        zip_buffer.seek(0)
        return zip_buffer.getvalue(), archive_name

if __name__ == "__main__":
    print("test")