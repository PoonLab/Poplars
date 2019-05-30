"""
Wrapper script to ensure that MAFFT use system-specific binaries
MAFFT versions: 7.427 (Linux), 7.482 (Windows), 7.427 (Mac OS X)
"""

import sys
import os
import subprocess
import logging
import argparse


def run_mafft(file_path):
    """
    Runs MAFFT on Linux, Windows, or Mac OS X
    :param file_path: path to the input file
    :return output: the output of MAFFT
    """

    try:
       sys.platform.startswith("linux") or sys.platform("win") or sys.platform == "darwin"

    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    # Path to find MAFFT
    script_path = os.path.dirname(os.path.abspath(__file__))

    if sys.platform.startswith("linux"):
        bin_path = os.path.join(script_path, 'bin/mafft-linux64/mafft.bat')

    elif sys.platform.startswith("win"):
        bin_path = os.path.join(script_path, 'bin/mafft-win/mafft.bat')

    else:
        bin_path = os.path.join(script_path, 'bin/mafft-mac/mafft.bat')

    if not os.path.isfile(bin_path):
        logging.error("No file exists.")

    raw_output = subprocess.check_output([bin_path, '--quiet', file_path])

    return raw_output


def parse_args():
    parser = argparse.ArgumentParser(
        description="Wrapper to run MAFFT"
    )
    parser.add_argument("infile",
                        help="Path to the input file")

    return parser.parse_args()


def main():
    args = parse_args()
    run_mafft(args.infile)


if __name__ == '__main__':
    main()
