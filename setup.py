#!/usr/bin/env python3

from distutils import log
from distutils.core import setup
from setuptools.command.install import install
import os


class OverrideInstall(install):
    """
    Adapted from:
    https://stackoverflow.com/questions/5932804/set-file-permissions-in-setup-py-file/25761434
    """
    def run(self):
        uid, gid = 0, 0  # root user
        mode = 0o755
        set_data_dir = False
        install.run(self)


mafft_files = []
for root, dirs, files in os.walk('bin/mafft-linux64/'):
    if len(files) > 0:
        mafft_files.append( (root, files) )

setup(
    name='poplars',
    packages=['poplars'],
    package_data={'poplars': [
        'bin/mafft-linux64/mafft.bat',
        'bin/mafft-linux64/mafftdir/bin/mafft',
        'bin/mafft-linux64/mafftdir/libexec/*'
    ]},
    cmdclass={'install': OverrideInstall}
)

