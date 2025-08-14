from __future__ import annotations

import argparse
import shlex
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import IO

from Bio.SeqUtils import MeltingTemp as mt
from pandas import DataFrame

from ..RNAProbesUtil import BufferedProgramObject, ProgramObject, run_command_line
from ..util import (path_string, validate_arg, parse_file_input,
                      DiscontinuousRange, input_range, validate_doesnt_throw, input_path, input_path_string, path_arg)
from ..RNAUtil import CT_to_sscount_df

undscr = ("->" * 40)
copyright_msg = ("\n" * 5) + (" \x1B[3m TFOFinder\x1B[0m  Copyright (C) 2025 Avi Kohn, 2022  Irina E. Catrina\n"
      "This program comes with ABSOLUTELY NO WARRANTY;\n"
      "This is free software, and you are welcome to redistribute it\n"
      "under certain conditions; for details please read the LICENSE.txt file.\n\n"
      "Feel free to use the CLI or to run the program directly with command line arguments \n"
      "(view available arguments with --help).\n\n"
      f"{undscr}\n\n"
      "WARNING: Previous files will be overwritten or appended!  Save them in a\n"
      "different location than the current input file, or rename them.\n\n"
      f"{undscr}")

probeMin = 4
probeMax = 30 #inclusive
exported_values = {"probeMin": probeMin, "probeMax": probeMax}


base_complement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})


def validate_arguments(probe_lengths: str, filename="", **ignore) -> dict:
    validate_arg(parse_file_input(filename).suffix == ".ct", "The given file must be a valid .ct file")
    probe_length_range = validate_doesnt_throw(DiscontinuousRange.template(probeMin, probeMax+1, True), probe_lengths, msg=f"The probe length is not between {probeMin} and {probeMax} inclusive")
    return {"probe_lengths": probe_length_range}

#todo: add support for sscount csv

def calculate_result(filein: IO[str], probe_lengths: DiscontinuousRange, filename: str="", arguments: argparse.Namespace =None, output_dir: Path = None) -> ProgramObject:
    """
    Find TFO probes
    This method should be inside a with expression in order to close the given file
    :param filein: the file object of the file to be read
    :param probe_lengths: the length of TFO probes to look for
    :param filename: the name of the inputted file
    :param arguments: a dictionary that controls printing during execution
    :return: str
    """

    fname = parse_file_input(filename).stem
    program_object = get_program_object(fname, arguments, output_dir)
    sscount_df, structure_count = CT_to_sscount_df(filein, True, program_object.save_buffer("[fname]_sscount.csv"))
    #todo: ask if can get rid of this and place before

    if should_print(arguments): print('Number of Structures = ' + str(structure_count) + ' \n\n...Please wait...\n')
    #temp
    all_probes_by_length = get_consecutive_not_ss(probe_lengths, sscount_df)
    program_object.set_result_args(final_result=get_final_string(filename, probe_lengths, structure_count, all_probes_by_length, sscount_df))

    # todo: ask if should add extra newline at end (trivial issue)
    with program_object.open_buffer(f"[fname]_TFO_probes.txt", "w" if arguments.overwrite else 'a') as buffer:
        buffer.write(program_object.result_obj.final_result)
    return program_object

def parse_arguments(args: str | list, from_command_line = True):
    args = get_argument_parser().parse_args(args if isinstance(args, list) else shlex.split(args))
    args.from_command_line = from_command_line  # denotes that this is from the command line
    return args

def run(args="", from_command_line = True):
    arguments = parse_arguments(args, from_command_line = from_command_line)
    if should_print(arguments):
       print(copyright_msg)

    filein = input_path_string(".ct", msg='Enter the ct file path and name: ', fail_message='This file is invalid (either does not exist or is not a .ct file). Please input a valid .ct file: ',
                        initial_value= arguments.file, retry_if_fail=arguments.from_command_line)
    mb_userpath, fname, suffix = parse_file_input(filein, arguments.output_dir)

    # probe_length = input_int_in_range(min=probeMin, max=probeMax + 1,
    #                                   msg=f"Enter the length of TFO probe; a number between {probeMin} and {probeMax} inclusive: ",
    #                                   fail_message=f'You must type a number between {probeMin} and {probeMax}, try again: ',
    #                                   initial_value=arguments.probe_length, retry_if_fail=arguments.from_command_line)
    probe_length = input_range(min=probeMin, max=probeMax + 1, force_increasing=True,
                                      msg=f"Enter the length of TFO probe; a number or range between {probeMin} and {probeMax} inclusive: ",
                                      fail_message=f'You must type a number (or a range) between {probeMin} and {probeMax}, try again: ',
                                      initial_value=arguments.probe_length, retry_if_fail=arguments.from_command_line)
    # convert the ct file to a txt
    with open(filein, 'r') as file:
        program_object = calculate_result(file, probe_length, filename=filein, arguments=arguments,
                                             output_dir=mb_userpath)

    if should_print(arguments, is_content_verbose=True):
        from textwrap import indent
        print("Results: \n" + indent(program_object.result_obj.final_result, " " * 4))
    if should_print(arguments):
        print("Calculation complete. Find the result in " + str(mb_userpath / f"{fname}_TFO_probes.txt"))

#region************************************************************************************************
#******************************************************************************************************
#**************************************   Util Functions   ********************************************
#******************************************************************************************************
#******************************************************************************************************

def get_program_object(file_stem, arguments, output_dir=None):
    constructor = ProgramObject if arguments.from_command_line else BufferedProgramObject
    return constructor(output_dir, file_stem, arguments)

def should_print(arguments, is_content_verbose = False):
    return arguments and arguments.from_command_line and not arguments.quiet and (not is_content_verbose or arguments.verbose)

def get_consecutive_not_ss(probe_lengths: DiscontinuousRange, sscount_df : DataFrame) -> Iterable[tuple[int, DataFrame]]:
    # get the double stranded elements with bse A or G
    dscount = sscount_df.loc[(sscount_df.base.isin(["A", "G"]) & (sscount_df.sscount != 20))]
    dscount = dscount.drop('sscount', axis=1)

    # get the first element of a 9 length probe
    consec = dscount.loc[dscount.baseno.diff(-1) == -1]  # and the one after
    return ((length, consec.loc[consec.baseno.diff(-(length - 2)) == -(length - 2)]) for length in reversed(probe_lengths))



def get_final_string(file_name : str, probe_lengths: DiscontinuousRange, structure_count: int, consec: Iterable[tuple[int, DataFrame]], sscount_df):
    to_return = 'Results for ' + file_name + ' using ' + str(probe_lengths) + ' as parallel TFO probe length(s)\n' + \
                'Start Position,%GA,sscount,Parallel TFO Probe Sequence,Tm,Probe Length\n'
    to_return += "\n".join((get_final_string_section(structure_count, df, length, sscount_df) for length, df in consec if len(df) != 0))
    return to_return + "\n"

def get_final_string_section(structure_count: int, consec: DataFrame, probe_length: int, sscount_df):
    get_probe_result = lambda baseno: ",".join(map(str, sequence_probe(baseno, probe_length, structure_count, sscount_df))) + f",{probe_length}"
    return "\n".join(map(get_probe_result, consec.baseno))
# def seqTarget(df: DataFrame): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
#     max_base = df.base.iat[-1]
#     seq = ''.join(df.base)
#     size = len(df)
#     return df.sscount, df.baseno, max_base, df.base, seq, size

def parallel_complement(seq : str, complement = base_complement): #generate RNA complement
    return seq.translate(complement)[::1]

def sequence_probe(baseno: int, probe_len: int, structure_count: int, sscount_df: DataFrame):
    """
    Sequence one specific probe
    :param baseno: an integer representing the start of the probe, 1-indexed
    :param probe_len: the length of the probe
    :param sscount_df: the sscount dataframe
    :return:
    """
    index = baseno - 1
    probe = "".join(sscount_df.base[index: index + probe_len])
    complement = parallel_complement(probe)

    tml = int(mt.Tm_NN(complement, dnac1=50000, dnac2=50000, Na=100, nn_table=mt.RNA_NN1, saltcorr=1))

    per = int((probe.count('A') + probe.count('G')) / probe_len * 100)

    probe_sscounts = sscount_df.sscount[index:index + probe_len]
    avg_sscount = probe_sscounts.sum() / (probe_len * structure_count)
    return (baseno, per, avg_sscount, complement, tml) #returns baseno so it can easily be written to file

argument_parser = None
def get_argument_parser():
    global argument_parser
    if argument_parser is None:
        argument_parser = create_arg_parser()
    return argument_parser

def create_arg_parser():
    import argparse, functools
    parser = argparse.ArgumentParser(
        prog='TDOFinder',
        description='Triplex-forming oligonucleotide (TFO) target designer for RNA.')
    parser.add_argument("-f", "--file", type=path_string)
    parser.add_argument("-o", "--output-dir", type=functools.partial(path_arg, suffix=""))
    parser.add_argument("-p", "--probe-length", type=DiscontinuousRange.template(probeMin, probeMax+1, True),
                        metavar=f"[{probeMin}-{probeMax}]",
                        help=f'The length range of TFO probes, {probeMin}-{probeMax} inclusive')
    parser.add_argument("-s", "--emit-sscount", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-w", "--overwrite", action="store_true",
                        help="Overwrite the returned TFO Probes file. Default is to append")
    return parser

#******************************************************************************************************
#******************************************************************************************************
#***********************************   End of Util Functions   ****************************************
#******************************************************************************************************
#endregion*********************************************************************************************

if __name__ == "__main__":
    #if we succeed somehow (throught pythonpath, etc)...
    run_command_line(run, sys.argv[1:])