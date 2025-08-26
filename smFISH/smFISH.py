from __future__ import annotations

import datetime
import os, sys
from argparse import Namespace
import shlex

import pandas as pd
import itertools
import math
from pathlib import Path
from pandas import DataFrame, Series

from ..RNAProbesUtil import run_command_line, ProgramObject
from ..RNAUtil import RNAStructureWrapper, get_ct_nucleotide_length
from ..smFISH.ReverseDijkstra import ReverseDijkstra
from ..util import path_string, path_arg, input_bool, validate_arg, parse_file_input, input_path_string, \
    format_timedelta, validate_doesnt_throw

undscr = ("->" * 40) + "\n"
copyright_msg = (("\n" * 6) +
          f'smFISH_HybEff program  Copyright (C) 2025 Avi Kohn, 2022  Irina E. Catrina\n' +
          'This program comes with ABSOLUTELY NO WARRANTY;\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; for details please read the LICENSE.txt file.\n\n' +
          "Feel free to use the CLI or to run the program directly with command line arguments \n" +
          "(view available arguments with --help).\n\n" +
          undscr +
          "\nWARNING: Previous files will be overwritten or appended!\n" +
          undscr)

#These are constants that users may want to modify in order to get different probes (so are lowercase):
probe_length = 20
min_hybeff = .6
count_weight, hybeff_weight = .5, 7 #change this to modify how the values are weighed. Rarely matters, but sometimes does
def hybeff_modifier(hybeff) -> float:
    modified_value = (hybeff - min_hybeff) / (1 - min_hybeff) #convert [min_hybeff,1] -> [0,1]
    return modified_value*modified_value*modified_value*modified_value

# count_weight, hybeff_weight = 1, 7 #change this to modify how the values are weighed. Rarely matters, but sometimes does
# def hybeff_modifier(hybeff) -> float:
#     return hybeff * hybeff * hybeff * hybeff * hybeff

# count_weight, hybeff_weight = 1, 1 #change this to modify how the values are weighed
# def hybeff_modifier(hybeff) -> float:
#     return hybeff


#region #Constants:
PROBE_RETURN_COUNT = 48
PRECISION = 10
CONCENTRATION = 0.25e-6
GAS_CONSTANT = 0.001987
TEMP_K = 310.15 #37 C or 98.6 F
# so there's a limit on what a webserver will allow
IS_WEBAPP = os.environ.get("IS_WEB_APP")
MAX_WEBAPP_NUC_LENGTH = 4 * 1000 if IS_WEBAPP else 50 * 1000 #if > 50k, will take ~10 hours


COLS_TO_SAVE = ('Pos', "Oligo(5'->3')", 'Overall (kcal/mol)', 'Tm-Dup (degC)', 'Hybeff', 'fGC')
#endregion

exported_values = dict(maxWebappLength=MAX_WEBAPP_NUC_LENGTH) #max file size: 2mb if web app. OligoWalk is O(n^3) and bifold is also bad,

def validate_arguments(file_path: Path, arguments: Namespace, **ignore) -> dict:
    validate_arg(parse_file_input(file_path).suffix == ".ct", "The given file must be a valid .ct file")
    validate_arg(Path(file_path).exists(), msg="The ct file must exist")
    nuc_length = validate_doesnt_throw(get_ct_nucleotide_length, file_path, msg="The given CT file is invalid. Can't read the CT file.")
    validate_arg(nuc_length < MAX_WEBAPP_NUC_LENGTH, f"The RNA length must be below {MAX_WEBAPP_NUC_LENGTH} nucleotides "
                                                                              f"{'when using a webapp. Feel free to run the program, downloaded through our GitHub repository, on your own system' if IS_WEBAPP else 'when running the program. Feel free to change it manually, but it may take incredibly long'}")
    validate_arg(hasattr(arguments, 'intermolecular') and arguments.intermolecular is not None, "You must use decide whether to use intermolecular or not")
    return dict(nucleotide_length=nuc_length)


def calculate_result(file_path : str | Path, arguments: Namespace, output_dir: Path = None, **ignore) -> ProgramObject:
    output_dir, fname, _ = parse_file_input(file_path, output_dir or arguments.output_dir)
    get_missing_arguments(arguments)
    program_object = ProgramObject(output_dir=output_dir, file_stem=fname, arguments=arguments)
    probes = get_best_possible_probe_set(file_path, program_object)
    best_48 = get_best_probes(probes, program_object, count=PROBE_RETURN_COUNT)
    try_intermolecular(best_48, program_object)

    return program_object

def parse_arguments(args: str | list, from_command_line = True) -> Namespace:
    args = get_argument_parser().parse_args(args if isinstance(args, list) else shlex.split(args))
    args.from_command_line = from_command_line  # denotes that this is from the command line
    args.intermolecular = args.force_intermolecular or args.hybeff or None
    return args

def run(args="", from_command_line = True):
    arguments = parse_arguments(args, from_command_line=from_command_line)
    if should_print(arguments): print(copyright_msg)

    ct_filein = input_path_string(".ct", msg='Enter the ct file path and name: ', fail_message='This file is invalid (either does not exist or is not a .ct file). Please input a valid .ct file: ',
                        initial_value= arguments.file, retry_if_fail=arguments.from_command_line)

    program_object = calculate_result(ct_filein, arguments)

    if arguments.intermolecular:
        #no filtered_file??
        print(program_object.format_relative_path(f"Check the [fname]_best_{PROBE_RETURN_COUNT}_probes.csv for proposed smFISH probes. However, if not enough probes have been"
              +" selected given the initial selection criteria or only the CDS is targeted, please review the [fname]_best_probes_set.csv and [fname]_possible_matching_probes.csv to "
              +"select additional probes. Moreover, the intermolecular interactions of the probes should be taken into acocunt. Please review the [fname]_combined_output.csv file, and eliminate any probes with "
              + "intermolecular hybdridization free energy change < -10kcal/mol."))
    else:
        print(program_object.format_relative_path(f"Check the [fname]_best_{PROBE_RETURN_COUNT}_probes.csv for proposed smFISH probes. However, if not enough probes have been "
              "selected given the initial selection criteria or only the CDS is targeted, please review the [fname]_best_probes_set.csv "
              "and [fname]_possible_matching_probes.csv to select additional probes."))

#region************************************************************************************************
#******************************************************************************************************
#**************************************   Util Functions   ********************************************
#******************************************************************************************************
#******************************************************************************************************

def get_missing_arguments(arguments: Namespace) -> None:
    arguments.intermolecular = input_bool(msg="Do you want to run smFISH in intermolecular mode? (y/n)",
                                          initial_value=arguments.intermolecular,
                                          retry_if_fail=arguments.from_command_line)

def try_intermolecular(df_filtered: DataFrame, program_object: ProgramObject):
    if program_object.arguments.intermolecular:
        oligos = df_filtered['Oligo(5\'->3\')'].tolist()

        # Process the second program using the extracted oligos
        process_oligos(oligos, program_object)


def parse_file_name(filein, output_dir:Path=None):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    output_dir = output_dir or filein_path.parent  # Use the parent directory of the input file to save all files
    return output_dir, filein_path.stem, filein_path.suffix

def count_c_g(cell):
    return cell.count('C') + cell.count('G')


def can_have_path(val, first_match, possible_node):
    to_return = (possible_node.Pos >= val.Pos + probe_length + 2 if val is not None else True) #can't be within 22 (probe_length + 2), unless this is the root
    if to_return and first_match is not None:
        to_return = to_return and first_match.Pos <= possible_node.Pos <= first_match.Pos + probe_length + 2 #if it doesn't match this, we can choose first_match AND this
    return to_return

def alg_cost_mapper(val: Series):
    return count_weight + hybeff_weight * hybeff_modifier(val.Hybeff)

def get_filtered_df(df: DataFrame, program_object: ProgramObject) -> DataFrame:
    alg = ReverseDijkstra(df, value_mapper= alg_cost_mapper,
                            can_have_path=can_have_path, indexer = lambda df, i: df.iloc[i])
    max_val, path = alg.run()
    filtered_df = DataFrame(path)
    filtered_df.reset_index(drop=True, inplace=True)

    return filtered_df

def get_best_possible_probe_set(filein, program_object: ProgramObject) -> DataFrame:
    if not program_object.arguments.csv_file:
        matching_probes = get_matching_probes(filein, program_object)
    else:
        matching_probes = pd.read_csv(program_object.arguments.csv_file)
        validate_arg(set(COLS_TO_SAVE).issubset(set(matching_probes.columns)), "The csv file is invalid. It must contain the column(s): " + ", ".join(set(COLS_TO_SAVE).difference(set(matching_probes.columns))))
        #could also verify the datatypes, but this much should be fine. Should just throw if invalid, which is OK

    filtered_df = get_filtered_df(matching_probes, program_object)
    filtered_df.to_csv(program_object.save_buffer(f"[fname]_best_probes_set.csv"), index=False, float_format=f'%.{PRECISION}g')

    return filtered_df

def get_best_probes(df: DataFrame, program_object: ProgramObject, count=PROBE_RETURN_COUNT) -> DataFrame:
    probes = df.sort_values(["Hybeff"], ignore_index=True, kind="stable").head(count)
    probes.to_csv(program_object.save_buffer(f"[fname]_best_{count}_probes.csv"), index=False, float_format=f'%.{PRECISION}g')
    if should_print(program_object.arguments): print(f"{len(df)} probes found." + (f" Choosing the best {count}" if len(df) > count else ""))
    return probes

def equilibrium_constant(input):
    return math.e ** (-(input / (GAS_CONSTANT*TEMP_K)))

def get_size_warning(length: int):
    if length > 4000:
        estimated_seconds = (10 * 60 / (12000 ** 3) * (length ** 3))
        return f"Calculation may take a while due to RNA length. Estimated time to completion: {format_timedelta(datetime.timedelta(seconds=estimated_seconds), include_seconds=False)}"
    return ""

def get_matching_probes(filein: str, program_object: ProgramObject):
    if should_print(program_object.arguments):
        length = get_ct_nucleotide_length(filein)
        print(get_size_warning(length))
    df = RNAStructureWrapper.oligowalk(Path(filein),
                                       arguments=f"--structure -d -l {probe_length} -c {CONCENTRATION} -m 1 -s 3 --no-header",
                                       path_mapper=program_object.file_path,
                                       remove_input=program_object.arguments.delete_ct)
    # todo: ummm, 0.1 * 10???
    dG1FA, dG2FA, dG3FA = (df['Duplex (kcal/mol)'] + 0.2597 * 10,
                           df['Intra-oligo (kcal/mol)'] + 0.1000 * 10,
                           df['Break-Target (kcal/mol)'] + (0.0117 * abs(df['Break-Target (kcal/mol)'])) * 10)
    Koverall = (equilibrium_constant(dG1FA) /
                ((1 + equilibrium_constant(dG2FA)) * (1 + equilibrium_constant(dG3FA))))
    k_overall = CONCENTRATION * Koverall
    df['Hybeff'] = k_overall / (1 + k_overall)

    df['fGC'] = (df['Oligo(5\'->3\')'].apply(
        count_c_g)) / probe_length  # Apply the function to each cell in the DataFrame; GC fraction in each sequence
    df.rename(columns={'Pos.': 'Pos'}, inplace=True)

    df_filtered = df[(df.fGC >= 0.45) & (df.fGC <= 0.60) & (
                df.Hybeff >= 0.6)]  # & (df2.Pos >= 434) & (df2.Pos <= 1841)] #only CDS for oskRC
    df_filtered.reset_index(drop=True, inplace=True)

    df_cols_removed = df_filtered[list(COLS_TO_SAVE)]
    df_cols_removed.to_csv(program_object.save_buffer("[fname]_possible_matching_probes.csv"), sep=',', index=None) #not using float_format so it can be used as an input
    return df_cols_removed

def process_oligos(oligos: list, program_object: ProgramObject):
    pairs = get_pairs_file(oligos, program_object.file_path("[fname]_pairs.txt"))

    #run bifold with a dummy file
    energy_values = RNAStructureWrapper.bifold(f"[fname]_pairs.txt", f"[fname]somefile", f"[fname]pairs.out", program_object.file_path,
                                               remove_input=True)

    # Combine the pairs with the energy values
    with program_object.open_buffer("[fname]_combined_output.csv", 'w') as f:
        f.write("Seq#1,Seq#2,DG\n")
        for i in range(len(energy_values)):
            f.write(f"{pairs[i][0]},{pairs[i][1]},{energy_values[i]}\n")

def get_pairs_file(oligos: list, path: Path) -> list:
    pairs = list(itertools.combinations(oligos, 2))  # Convert to list for indexing
    with open(path, "w") as f:
        for a, b in pairs:
            f.write(a + ' ' + b + '\n')
    return pairs

argument_parser = None
def get_argument_parser():
    global argument_parser
    if argument_parser is None:
        argument_parser = create_arg_parser()
    return argument_parser
def create_arg_parser():
    import argparse, functools
    parser = argparse.ArgumentParser(
        prog='smFISH',
        description='Probe design for single-molecule fluorescence in situ hybridization (smFISH), considering secondary structures.')

    arg_group = parser.add_argument_group('Input file',
                                          'The input file to use')
    group = arg_group.add_mutually_exclusive_group()
    group.add_argument("-f", "--file", type=functools.partial(path_string, suffix=".ct"), help="The default file input.")
    group.add_argument("-cf", "--csv-file", type=functools.partial(path_string, suffix=".csv"),
                       help='For speed if smFISH was previously ran on a ct file. Input a previously outputted "*_possible_matching_probes.csv"')  # for testing, speed

    parser.add_argument("-o", "--output-dir", type=functools.partial(path_arg, suffix=""))
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-d", "--delete-ct", action="store_true", help="Remove the ct input file. Not recommended unless running from a server")

    arg_group = parser.add_argument_group('Intermolecular',
                                          'Intermolecular command line settings. If none given, will ask')
    group = arg_group.add_mutually_exclusive_group()
    group.add_argument("-i", "--intermolecular", dest="force_intermolecular", action="store_const", const="y")
    group.add_argument("-ni", "--no-intermolecular", "--hybeff", dest = "hybeff", action="store_const", const="n")

    return parser

def should_print(arguments: Namespace, is_content_verbose = False):
    return arguments and arguments.from_command_line and not arguments.quiet and (not is_content_verbose or arguments.verbose)

#******************************************************************************************************
#******************************************************************************************************
#***********************************   End of Util Functions   ****************************************
#******************************************************************************************************
#endregion*********************************************************************************************

if __name__ == "__main__":
    #if we succeed somehow (throught pythonpath, etc)...
    run_command_line(run, sys.argv[1:])