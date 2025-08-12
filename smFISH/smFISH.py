from __future__ import annotations

import os, sys
from argparse import Namespace
import shlex

import pandas as pd
import itertools
import math
from pathlib import Path
from pandas import DataFrame, Series

from ..RNASuiteUtil import run_command_line, ProgramObject
from ..RNAUtil import RNAStructureWrapper
from ..smFISH.ReverseDijkstra import ReverseDijkstra
from ..util import path_string, path_arg, input_bool, validate_arg, parse_file_input, input_path_string

undscr = ("->" * 40) + "\n"
copyright_msg = (("\n" * 6) +
          f'smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina\n' +
          'This program comes with ABSOLUTELY NO WARRANTY;\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; for details please read the GNU_GPL.txt file.\n\n' +
          "Feel free to use the CLI or to run the program directly with command line arguments \n" +
          "(view available arguments with --help).\n\n" +
          undscr +
          "\nWARNING: Previous files will be overwritten or appended!\n" +
          undscr)
# so there's a limit on what a webserver will allow
probe_length = 20
concentration = 0.25e-6
gas_constant = 0.001987
temp_K = 310.15 #37 C or 98.6 F
min_Hybeff = .6
max_webapp_size = 2*1000*1000 if os.environ.get("IS_WEB_APP") else 20*1000*1000 #alg is n^3, so is 1000 larger so basically infinite (if larger, won't even run)

count_weight, hybeff_weight = .5, 7 #change this to modify how the values are weighed. Rarely matters, but sometimes does
def hybeff_modifier(hybeff) -> float:
    modified_value = (hybeff - min_Hybeff) / (1 - min_Hybeff) #convert [min_Hybeff,1] -> [0,1]
    return modified_value*modified_value*modified_value*modified_value

exported_values = dict(maxWebappSize=max_webapp_size) #max file size: 2mb if web app. OligoWalk is O(n^3) and bifold is also bad,

# count_weight, hybeff_weight = 1, 7 #change this to modify how the values are weighed. Rarely matters, but sometimes does
# def hybeff_modifier(hybeff) -> float:
#     return hybeff * hybeff * hybeff * hybeff * hybeff

# count_weight, hybeff_weight = 1, 1 #change this to modify how the values are weighed
# def hybeff_modifier(hybeff) -> float:
#     return hybeff

def validate_arguments(file_path: Path, arguments: Namespace, **ignore) -> dict:
    validate_arg(parse_file_input(file_path).suffix == ".ct", "The given file must be a valid .ct file")
    validate_arg(Path(file_path).exists(), msg="The ct file must exist")
    validate_arg(os.path.getsize(file_path) < max_webapp_size, f"The CT file must be below {max_webapp_size / 1000 / 1000} MB when using a webapp")
    validate_arg(hasattr(arguments, 'intermolecular') and arguments.intermolecular is not None, "You must use decide whether to use intermolecular or not")
    return dict()


def calculate_result(file_path : str | Path, arguments: Namespace, output_dir: Path = None, **ignore):
    output_dir, fname, _ = parse_file_input(file_path, output_dir or arguments.output_dir)
    get_missing_arguments(arguments)
    program_object = ProgramObject(output_dir=output_dir, file_stem=fname, arguments=arguments)
    df_filtered = process_ct_file(file_path, program_object)

    try_intermolecular(program_object, df_filtered)

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

    calculate_result(ct_filein, arguments)

    if arguments.intermolecular:
        #no filtered_file??
        print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been"
              +" selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv and *_unfiltered_probes.csv to "
              +"select additional probes. Moreover, the intermolecular interactions of the probes should be taken into acocunt. Please review the *combined_output.csv file, and eliminate any probes with "
              + "intermolecular hybdridization free energy change < -10kcal/mol.")
    else:
        print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been "
              "selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv "
              "and *_unfiltered_probes.csv to select additional probes.")

#region************************************************************************************************
#******************************************************************************************************
#**************************************   Util Functions   ********************************************
#******************************************************************************************************
#******************************************************************************************************

def get_missing_arguments(arguments: Namespace) -> None:
    arguments.intermolecular = input_bool(msg="Do you want to run smFISH in intermolecular mode? (y/n)",
                                          initial_value=arguments.intermolecular,
                                          retry_if_fail=arguments.from_command_line)

def try_intermolecular(program_object: ProgramObject, df_filtered):
    if program_object.arguments.intermolecular:
        oligos = df_filtered['Oligo(5\'->3\')'].tolist()

        # Process the second program using the extracted oligos
        process_list_file(program_object, oligos)


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

def get_filtered_file(df: DataFrame, program_object: ProgramObject) -> DataFrame:
    alg = ReverseDijkstra(df, value_mapper= alg_cost_mapper,
                            can_have_path=can_have_path, indexer = lambda df, i: df.iloc[i])
    max_val, path = alg.run()
    filtered_df = DataFrame(path)
    filtered_df.reset_index(drop=True, inplace=True)

    #filtered_df.to_csv(program_object.save_buffer("[fname]filtered_file.csv"), index=False)  # Save the result to filtered_file.csv
    return filtered_df

def process_ct_file(filein, program_object: ProgramObject) -> DataFrame:
    arguments = program_object.arguments
    if not program_object.arguments.csv_file:
        matching_probes = get_matching_probes(filein, program_object)
    else:
        matching_probes = pd.read_csv(program_object.arguments.csv_file)

    filtered_df = get_filtered_file(matching_probes, program_object) #idk why only if not intermolecular, but whatever
    filtered_df.to_csv(program_object.save_buffer(f"[fname]_final_filtered_file.csv"), index=False, float_format='%.10g')

    return filtered_df

def equilibrium_constant(input):
    return math.e ** (-(input / (gas_constant*temp_K)))

def get_matching_probes(filein: str, program_object: ProgramObject):
    df = RNAStructureWrapper.oligowalk(Path(filein),
                                       arguments=f"--structure -d -l {probe_length} -c {concentration} -m 1 -s 3 --no-header",
                                       path_mapper=program_object.file_path,
                                       remove_input=program_object.arguments.delete_ct)
    # todo: ummm, 0.1 * 10???
    dG1FA, dG2FA, dG3FA = (df['Duplex (kcal/mol)'] + 0.2597 * 10,
                           df['Intra-oligo (kcal/mol)'] + 0.1000 * 10,
                           df['Break-Target (kcal/mol)'] + (0.0117 * abs(df['Break-Target (kcal/mol)'])) * 10)
    Koverall = (equilibrium_constant(dG1FA) /
                ((1 + equilibrium_constant(dG2FA)) * (1 + equilibrium_constant(dG3FA))))
    k_overall = concentration * Koverall
    df['Hybeff'] = k_overall / (1 + k_overall)

    df['fGC'] = (df['Oligo(5\'->3\')'].apply(
        count_c_g)) / probe_length  # Apply the function to each cell in the DataFrame; GC fraction in each sequence
    df.rename(columns={'Pos.': 'Pos'}, inplace=True)

    df_filtered = df[(df.fGC >= 0.45) & (df.fGC <= 0.60) & (
                df.Hybeff >= 0.6)]  # & (df2.Pos >= 434) & (df2.Pos <= 1841)] #only CDS for oskRC
    df_filtered.reset_index(drop=True, inplace=True)

    df_cols_removed = df_filtered[['Pos', "Oligo(5'->3')", 'Overall (kcal/mol)', 'Tm-Dup (degC)', 'Hybeff', 'fGC']]
    df_cols_removed.to_csv(program_object.save_buffer("[fname]_possible_matching_probes.csv"), sep=',', index=None,
                           float_format='%.10g')
    return df_cols_removed
def process_list_file(program_object: ProgramObject, oligos):
    output_dir = program_object.output_dir #temp
    fname = program_object.file_stem
    pairs = list(itertools.combinations(oligos, 2))  # Convert to list for indexing
    
    with open(output_dir / f"{fname}pairs.txt", 'w') as f:
        for a, b in pairs:
            f.write(a + ' ' + b + '\n')

    #run bifold with a dummy file
    energy_values = RNAStructureWrapper.bifold(f"[fname]pairs.txt", f"[fname]somefile", f"[fname]pairs.out", program_object.file_path,
                                               remove_input=True)

    # Combine the pairs with the energy values
    with program_object.open_buffer("[fname]combined_output.csv", 'w') as f:
        f.write("Seq#1,Seq#2,DG\n")
        for i in range(len(energy_values)):
            f.write(f"{pairs[i][0]},{pairs[i][1]},{energy_values[i]}\n")

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
    parser.add_argument("-f", "--file", type=functools.partial(path_string, suffix=".ct"))
    parser.add_argument("-cf", "--csv-file", type=functools.partial(path_string, suffix=".csv"), help="For testing. Input a previously outputted all probes csv file") #for testing, speed
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