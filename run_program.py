from __future__ import annotations
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from rnaprobes.PinMol import pinmol
from rnaprobes.RNASuiteUtil import run_command_line
from rnaprobes.TFOFinder import tfofinder
from rnaprobes.smFISH import smFISH
from rnaprobes.util import input_value


programs = {
    "tfofinder": tfofinder.run,
    "pinmol": pinmol.run,
    "smfish": smFISH.run
}
if __name__ == "__main__":
    program = input_value("Input a program (either tfofinder, pinmol, or smfish): ", str.lower,
                            lambda program: program in programs.keys(), retry_if_fail=len(sys.argv) >= 2,
                          initial_value=sys.argv[1].lower() if len(sys.argv) >= 2 else None)
    run = programs[program]
    run_command_line(run, sys.argv[2:])