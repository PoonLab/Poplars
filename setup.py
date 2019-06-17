#!/usr/bin/env python3
from distutils import log
from distutils.core import setup
from setuptools.command.install import install
import os
import sys
from setuptools.command.install_scripts import install_scripts

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
        for filepath in self.get_outputs():
            path = os.path.dirname(filepath)
            if path.endswith('ref_genomes'):
                log.info('Changing permissions of {0} to {1:o}'.format(filepath, mode))
                os.chmod(filepath, mode)


data_files = ['ref_genomes/K03455.fasta', 'ref_genomes/M33262.fasta', 
    'ref_genomes/HIV1_Mgroup.fasta']

# identify platform
if sys.platform.startswith("linux"):
    mafft_dir = 'bin/mafft-linux64/'
elif sys.platform.startswith("win"):
    mafft_dir = 'bin/mafft-win/'
elif sys.platform == "darwin":
    mafft_dir = 'bin/mafft-mac/'
else:
    sys.stderr.write("Warning: failed to recognize platform.")
    sys.exit()

# generate listing of MAFFT-associated files for user's platform
mafft_files = []
for root, dirs, files in os.walk('bin/mafft-linux64/'):
    if len(files) > 0:
        mafft_files.append((root, files))

setup(
    name='poplars',
    packages=['poplars'],
    package_data={'poplars': mafft_files + data_files},
    cmdclass={'install': OverrideInstall}
)


