from __future__ import annotations
import sys
from pathlib import Path

#sys.path.append(str(Path(__file__).resolve().parent.parent))

from .PinMol import pinmol
from .RNASuiteUtil import run_command_line
from .TFOFinder import tfofinder
from .smFISH import smFISH
from .util import input_value


programs = {
    "tfofinder": tfofinder.run,
    "pinmol": pinmol.run,
    "smfish": smFISH.run
}
def run(args: list):
    program = input_value("Input a program (either tfofinder, pinmol, or smfish): ", str.lower,
                          lambda program: program in programs.keys(), retry_if_fail=len(args) >= 1,
                          initial_value=args[0].lower() if len(args) >= 1 else None)
    run = programs[program]
    run_command_line(run, args[1:])
if __name__ == "__main__":
    run(sys.argv[1:])