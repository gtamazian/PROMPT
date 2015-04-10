#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
from os.path import join


class MinActionPathPDB(object):
    """
    Implements routines to process a PDB file produced by
    MinActionPath.
    """

    def __init__(self, filename):
        """
        Create an object from the specified PDB file.

        :param filename: a name of a PDB file produced by MinActionPath
        :type filename: str
        """
        self.__filename = filename
        self.models = []
        self.__read_pdb()

    def __read_pdb(self):
        """
        Read lines from the input PDB file and form a list of them
        which entries correspond to single models.

        :return: a list of lines from the PDB file, one element for
            one model
        :rtype: list
        """
        current_model = []
        with open(self.__filename) as pdb_file:
            for line in pdb_file:
                current_model.append(line)
                if line.startswith('ENDMDL'):
                    self.models.append(current_model)
                    current_model = []

    def reduce(self, m):
        """
        Reduce the number of intermediate configurations in the
        transformation to the specified value. The transformations
        are chosen to be uniformly distributes between the first and
        last one.

        :param m: the required number of intermediate configurations
        :type m: int
        """
        n = len(self.models)
        a = float(n - 1)/(m - 1)
        b = 1 - a
        new_models = []
        for i in xrange(m):
            new_models.append(self.models[int(a*(i + 1) + b) - 1])
        self.models = new_models

    def write_pdb_files(self, directory='.', prefix='pdb_'):
        """
        Write single-model PDB files to the specified directory.

        :param directory: a directory to write the PDB files to
        :param prefix: a prefix of PDB file names
        :type directory: str
        :type prefix: str
        """
        for i, model in enumerate(self.models):
            output_filename = join(directory, prefix)
            with open('{}{}.pdb'.format(output_filename, i), 'w') as \
                    output:
                for line in self.models[i]:
                    output.write(line)


def parse_command_line_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser('Split a PDB file produced by '
                                     'MinActionPath to single-model '
                                     'PDB files')
    parser.add_argument('pdb_file', help='a PDB file by MinActionPath')
    parser.add_argument('-d', '--directory', default='.',
                        help='a directory to write single-model PDB '
                             'files to')
    parser.add_argument('-p', '--prefix', default='pdb_',
                        help='a prefix of names of single-model PDB '
                             'files')
    parser.add_argument('-m', '--models', type=int,
                        help='the number of models which PDB files '
                             'will be output')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_command_line_args()
    transformation = MinActionPathPDB(args.pdb_file)
    if args.models:
        transformation.reduce(args.models)
    transformation.write_pdb_files(args.directory, args.prefix)
