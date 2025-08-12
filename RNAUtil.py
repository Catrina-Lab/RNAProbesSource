from __future__ import annotations

import json
import os
import platform
import shlex
import subprocess
from argparse import ArgumentError
from collections.abc import Callable
from os import PathLike
from platform import architecture

import pandas as pd
from pathlib import Path
from pandas import DataFrame
from typing import IO

from pandas._typing import WriteBuffer

from .util import remove_files, ValidationError

match = ["ENERGY", "dG"] #find header rows in ct file

def CT_to_sscount_df(file: IO[str], save_to_file: bool = None, output_file: Path = None) -> tuple[DataFrame, int]:
    ct_df, structure_count = convert_ct_to_dataframe(file)
    sscount_df = getSSCountDF(ct_df, save_to_file, output_file)
    return sscount_df, structure_count

def convert_ct_to_dataframe(file: IO[str]) -> tuple[DataFrame, int]:
    """
        Convert the ct file to a dataframe. This also gives the number of structures in the ct file.
        This method should be inside of a with expression.
        :param filename: the file to read and convert into a dataframe.
        :return: (Dataframe, int structure_count)
        """

    # dtype={"baseno": int, "base": str, "bs_bind": int}
    try:
        file.seek(0)
        ct_df = pd.read_csv(file, sep='\\s+', usecols=[0,1,4], names=["baseno", "base", "bs_bind"], engine='python', skiprows=1)
        initial_row_count = len(ct_df)
        ct_df = ct_df[~ct_df["base"].isin(match)]
        structure_count = initial_row_count - len(ct_df) + 1 #need to add 1 because we skip the first row (in order to get the correct column count)
        ct_df = ct_df.astype({"baseno": int, "base": str, "bs_bind": int})
        return ct_df, structure_count
    except Exception as e:
        raise ValidationError("Can't parse the CT file. Is the CT file invalid?") from e


def getSSCountDF(ct_dataframe : DataFrame, save_to_file: bool = None, output_file: str | PathLike[str] | WriteBuffer[bytes] | WriteBuffer[str] = None) -> DataFrame:
    """
    Get the SSCount from a dataframe
    :param ct_dataframe: ct converted to a dataframe
    :return: Dataframe
    """

    df_grouped = ct_dataframe.groupby(['baseno', 'base'], as_index = False).agg(lambda x: x[x == 0].count())
    df_grouped.rename(columns={'bs_bind': 'sscount'}, inplace=True)
    df_grouped = df_grouped.reindex(columns=['baseno','sscount','base'])


    df_grouped['base'] = df_grouped.base.replace('T', 'U')

    if save_to_file:
        df_grouped.to_csv(output_file, index=False, header=False)
    return df_grouped

def _map_all(path_mapper: Callable[[str], Path | str], *files: str | Path) -> tuple[Path | str, ...]:
    return tuple((path_mapper(file) if isinstance(file, str) else file) for file in files)

base_folder = Path(__file__).parent / "RNAStructure_Binaries"
folders = None
rna_structure_directory = None
def get_folders():
    global folders
    if folders is None:
        path = base_folder / "architectures.json"
        with open(path, "r") as file:
            folders = json.load(file)
    return folders
def get_RNAStructure_directory() -> Path:
    """
    Get the RNAStructure directory for optimized RNAStructure programs
    :return:
    """
    global rna_structure_directory
    rna_structure_directory = (rna_structure_directory or
                               base_folder / get_folders().get(f"{platform.system()},{platform.machine()}", "")) #use root if none found, will just be the system
    if rna_structure_directory != base_folder: os.environ["DATAPATH"] = str(base_folder / "data_tables")
    return rna_structure_directory

def get_program(program: str):
    directory = get_RNAStructure_directory()
    file = next((file for file in directory.iterdir() if file.stem.lower() == program.lower()), None)
    return file or program

#specific to RNAStructure
def _run_program(program, files_in: list | str | Path, file_out: str | Path, arguments: str = "", remove_input: bool = False):
    files_in = files_in if isinstance(files_in, list) else [files_in]
    try:
        program_path = get_program(program)
        subprocess.check_output([program_path, *files_in, file_out, *shlex.split(arguments)])
        return file_out
    except subprocess.CalledProcessError as e:
        print("SUBPROCESS ERROR:")
        print(e.stderr)
        raise Exception(f"Error when running program {program}: " + str(e.stderr))
    finally:
        if remove_input: remove_files(*files_in)

class RNAStructureWrapper:
    @staticmethod
    def oligoscreen(input: pd.Series, file_name: str, path_mapper: Callable[[str], Path | str] = lambda x: x, arguments: str = "") -> DataFrame:
        input_path = path_mapper(f"{file_name}_oligoscreen_input.lis")
        output_path = path_mapper(f"{file_name}_oligoscreen_output.csv")
        try:
            input.to_csv(input_path, index=False, header=False)

            _run_program("oligoscreen",  input_path, output_path, *shlex.split(arguments))
            read_oligosc = pd.read_csv(output_path, delimiter='\t', usecols=[1, 2, 3])
            return read_oligosc
        finally:
            remove_files(input_path, output_path)

    @staticmethod
    def fold(file_in: str | PathLike[str], file_out: str | PathLike[str], path_mapper: Callable[[str], Path | str] = lambda x: x,
                    arguments: str = "", remove_input: bool = False) -> Path | str:
        seq_file, ct_file = _map_all(path_mapper, file_in, file_out)
        _run_program("fold", seq_file, ct_file, arguments, remove_input=remove_input)
        return ct_file

    @staticmethod
    def draw(file_in: str | PathLike[str], file_out: str | PathLike[str], path_mapper: Callable[[str], Path | str] = lambda x: x,
                    arguments: str = "", remove_input: bool = False) -> Path | str:
        ct_file, svg_file = _map_all(path_mapper, file_in, file_out)
        _run_program("draw", ct_file, svg_file, arguments, remove_input=remove_input)

        return svg_file

    @staticmethod
    def oligowalk(file_in: str | Path, path_mapper: Callable[[str], Path] = lambda x: x,
                  arguments: str = "--structure -d -l 20 -c 0.25uM -m 1 -s 3", remove_input: bool = False,
                  keep_output = False) -> DataFrame:
        file_in, = _map_all(path_mapper, file_in)
        file_out = file_in.parent / (file_in.stem + "_oligowalk_output.txt")
        try:
            _run_program("OligoWalk", file_in, file_out, arguments, remove_input=remove_input)
            return pd.read_csv(file_out, skiprows=3, sep='\t')
        finally:
            if not keep_output: remove_files(file_out)

    @staticmethod
    def bifold(file_in1: str | Path, file_in2: str | Path, file_out: str | Path, path_mapper: Callable[[str], Path | str] = lambda x: x,
               arguments: str = "--DNA --intramolecular --list", remove_input: bool = False) -> list[str]:
        file_in1, file_in2, file_out = _map_all(path_mapper, file_in1, file_in2, file_out)
        try:
            _run_program("bifold", [file_in1, file_in2], file_out, arguments, remove_input=remove_input)
            with open(file_out, 'r') as f: energy_values = [line.strip() for line in f.readlines()]
            return energy_values
        finally:
            remove_files(file_out)
