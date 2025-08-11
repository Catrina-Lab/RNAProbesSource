from __future__ import annotations
from pathlib import Path
from unittest import TestCase
from pinmol import run
import pinmol
import shlex

class Test(TestCase):
    def test_program(self):
        output_file = Path(__file__).parent.parent.parent / "Temp_files" / "example_large.ct"
        run(shlex.split(fr'-p 20 -f "{output_file}" --start 1 --end -1 -w -nb'))
        path = Path(output_file)
        output_dir = path.parent
        verify_program_data(self, path.stem, output_dir=output_dir, reference_dir=output_dir / "pinmol_files" / "expected_pinmol_large_all_stable",
                            svg_files_output_dir=pinmol.svg_dir_name, svg_files_reference_dir="")

    def test_blast_program(self):
        output_dir = Path(__file__).parent.parent.parent / "Temp_files"
        output_file = output_dir / "example_large.ct"
        run(shlex.split(fr'-p 20 -f "{output_file}" --start 1 --end -1 -w -n -bf "{output_dir / "example_large_Alignment.xml"}"'))
        path = Path(output_file)
        verify_program_data(self, path.stem, output_dir=output_dir, reference_dir=output_dir / "pinmol_files" / "example_large_blast_stable",
                            svg_files_output_dir=pinmol.svg_dir_name, svg_files_reference_dir="")

def verify_program_data(tester: TestCase, file_stem: str, output_dir: Path, reference_dir: Path, svg_max = 50, svg_files_output_dir: str = "", svg_files_reference_dir: str = None):
    svg_files_reference_dir = svg_files_reference_dir or svg_files_output_dir
    regular_files = ["_blast_picks.fasta", "_DG_probes.csv", "_GC_probes.csv",
                     "_probes_sortedby5.csv", "_sscount.csv"] #check  "_Final_molecular_beacons.csv" manually
    file_names = [file_stem + file_name for file_name in regular_files]
    for file in file_names:
        assert_files_equal(tester, output_dir / file, reference_dir / file)
    for svg_num in range (1, svg_max+1):
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

def format_lines(lines, path, ignore_header_case = True):
    lines[-1] = lines[-1].rstrip()
    if path.suffix == ".csv" and ignore_header_case:
        lines[0] = lines[0].lower()
    return lines