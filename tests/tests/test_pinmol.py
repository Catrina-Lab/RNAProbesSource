from __future__ import annotations

import uuid
from pathlib import Path
from unittest import TestCase

import pytest

from ...PinMol.pinmol import run
from ...PinMol import pinmol
import shlex

from rnaprobes.util import safe_remove_tree

FILES = ("[fname]_Final_molecular_beacons.txt", "[fname]_best_probes.csv", "[fname]_blast_picks.fasta",
         "[fname]_all_probes_sortedby5.csv", "[fname]_GC_bounded_probes.csv", "[fname]_sscount.csv")
RUN_BLAST_FILES = ("[fname]_blast_result.xml",)

test_dir = Path(__file__).parent.parent
test_output_path = test_dir / "output" / "PinMol"
example_file_path = test_dir / "test_example_files"
test_file_path = test_dir / "test_example_files" / "PinMol"

class Test(TestCase):
    def test_large(self):
        run_test(self, "example_large", "no_blast/large",
                 r'-p 20 -f [fpath] --start 1 --end -1 -w -nb')
    def test_small(self):
        run_test(self, "example_small", "no_blast/small",
                 r'-p 20 -f [fpath] --start 1 --end -1 -w -nb')
    def test_super_large(self):
        run_test(self, "example_super_large", "no_blast/super_large",
                 r'-p 20 -f [fpath] --start 1 --end -1 -w -nb')

    def test_blast_program(self):
        run_test(self, "example_large", "blast/large",
                 fr'-p 20 -f [fpath] --start 1 --end -1 -w -bf "{test_file_path / "blast" / "example_large_blast_result.xml"}"')

    # def test_
class TestSlow(TestCase):
    @pytest.mark.slow
    def test_slow(self):
        run_test(self, "example_super_large", "no_blast/super_large",
                 r'-p 20 -f [fpath] --start 1 --end -1 -w -nb')


def run_test(tester: TestCase, file_stem: str, reference_dir_name: str,
             arguments: str, is_blast_run=False,
             example_path: Path = example_file_path):
    output_dir = test_output_path / str(uuid.uuid4())
    output_dir.mkdir(parents=True)
    try:
        run(shlex.split(arguments.replace("[fpath]", f'"{example_path / (file_stem + ".ct")}"') +
                        fr' -o "{str(output_dir.resolve())}"'), from_command_line=False)
        verify_program_data(tester, file_stem, output_dir=output_dir,
                            reference_dir=reference_dir_name, is_blast_run=is_blast_run)
    finally:
        safe_remove_tree(output_dir, test_output_path)


def verify_program_data(tester: TestCase, file_stem: str, output_dir: Path, reference_dir: str,
                        is_blast_run = False,
                        svg_max = 50, svg_files_output_dir: str = pinmol.svg_dir_name, svg_files_reference_dir: str = None):
    svg_files_reference_dir = svg_files_reference_dir or svg_files_output_dir
    reference_dir = test_file_path / reference_dir

    files = FILES + (RUN_BLAST_FILES if is_blast_run else tuple())
    files = [string.replace("[fname]", file_stem) for string in files]

    for file in files:
        assert_files_equal(tester, output_dir / file, reference_dir / file)

    for svg_num in range(1, svg_max+1):
        ref_svg = reference_dir / svg_files_reference_dir / f"{file_stem}_{svg_num}.svg"
        if ref_svg.exists():
            assert_files_equal(tester, output_dir / svg_files_output_dir / f"{file_stem}_{svg_num}.svg", ref_svg)

def assert_files_equal(tester: TestCase, test_path: Path, ref_path: Path):
    try:
        with open(test_path) as tst_file, open(ref_path) as ref_file:
            lines1 = format_lines(list(tst_file), test_path)
            lines2 = format_lines(list(ref_file), ref_path)
            tester.assertListEqual(
                lines1,
                lines2)
    except AssertionError as e:
        raise AssertionError(f"File {test_path} is different from the reference file {ref_path}") from e

def format_lines(lines, path, ignore_header_case = False):
    lines[-1] = lines[-1].rstrip()
    if path.suffix == ".csv" and ignore_header_case:
        lines[0] = lines[0].lower()
    return lines