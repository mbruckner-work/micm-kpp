"""
Copyright (C) 2024
National Center for Atmospheric Research,
SPDX-License-Identifier: Apache-2.0
(Software Package Data Exchange)

File:
    mozart_to_micm.py

Usage:
    python mozart_to_micm.py
    python mozart_to_micm.py --help
"""

import os
import sys
import argparse
import logging
import json
import yaml
from glob import glob

from parse_kpp_utils import is_float, parse_term


__version__ = 'v1.00'


unknowns = list()


def read_mozart_config(config_dir, config_name):
    """
    Read all Mozart config files in a directory

    Parameters
        (str) mozart_dir: Mozart directory

    Returns
        (list of str): all lines from all config files
    """

    suffixes = ['.spc', '.eqn', '.def']

    lines = list()

    for suffix in suffixes:
        files = glob(os.path.join(config_dir, config_name + '*' + suffix))
        logging.debug(files)
        for filename in files:
            with open(filename, 'r') as f:
                lines.extend(f.readlines())

    # remove empty lines and tabs
    lines = [line.replace('\t', '') for line in lines if line.strip()] 

    for line in lines:
        logging.debug(line.strip())

    return lines


def split_by_section(lines):
    """
    Split config lines by section

    Parameters
        (list of str) lines: all lines config files

    Returns
        (dict of list of str): lines in each section
    """

    sections = {'#Variable': [],
                '#Fixed': [],
                '#Equations': []}

    joined_lines = ''.join(lines)
    section_blocks = joined_lines.split('#')

    for section in sections:
        for section_block in section_blocks:
            if section.replace('#', '') in section_block:
                sections[section].extend(section_block.split('\n')[1:-1])

    return sections


def parse_species(lines, fixed=False, tolerance=1.0e-12):
    """
    Generate species JSON

    Parameters
        (list of str) lines: lines of species section
        (bool) fixed: set constant tracer
        (float) tolerance: absolute tolerance

    Returns
        (list of dict): list of MICM species entries
    """

    species_json = list() # list of dict

    for line in lines:
        if '->' in line:
            lhs, rhs = tuple(line.split('->'))
            logging.debug((lhs, rhs))
        else:
            lhs, rhs = line, None
            logging.debug(lhs)
        species_dict = {'name': lhs.replace('{', '').replace('}', '').strip().lstrip(),
            'type': 'CHEM_SPEC'}
        if rhs is not None:
            species_dict['__formula'] = rhs.replace('{', '').replace('}', '').strip().lstrip()
        if fixed:
            species_dict['tracer type'] = 'CONSTANT'
        else:
            species_dict['absolute tolerance'] = tolerance
        species_json.append(species_dict)

    return species_json


def parse_equations(lines):
    """
    Generate MICM equation JSON

    Parameters
        (list of str) lines: lines of equation section

    Returns
        (list of dict): list of MICM equation entries
    """

    equations = list() # list of dict

    joined_lines = '____'.join(lines) + '____'
    equation_sets = joined_lines.split('[')

    for equation_set in equation_sets:
        if ']' in equation_set:
            label, lines____ = tuple(equation_set.split(']'))
            label = label.strip().lstrip()
            logging.info('equation set ' + label)
            lines = lines____.split('____')[0:-1]
            parse_result = parse_equation_set(label, lines)
            if len(parse_result) > 0:
                equations.append(parse_result)

    return equations


def parse_equation_set(label, lines):
    """
    Generate MICM equation JSON

    Parameters
        (list of str) lines: lines of equation section

    Returns
        (list of dict): list of MICM equation entries
    """

    logging.debug('parse_equation_set ' + label)
    logging.debug('line count:' + str(len(lines)))

    equations = list() # list of dict

    line_sets = list()

    idx = -1
    for line in lines:
        logging.debug(line)
        if '->' in line:
            line_sets.append(list())
            idx += 1
        line_sets[idx].append(line)

    # print(line_sets)

    eqn_number = 1

    for lines in line_sets:
        rhs_combo = ''
        coeffs = ''

        for line in lines:
            logging.debug(line)

            # split on reaction delimiter into left hand and right hand sides 
            if '->' in line:
                lhs, rhs = tuple(line.split('->'))
                lhs, rhs = lhs.strip().lstrip(), rhs.strip().lstrip()
                if ';' in rhs:
                    rhs, coeffs = tuple(rhs.split(';'))
                    rhs = rhs.strip().lstrip()
                    coeffs = coeffs.replace(' ', '')
                else:
                    coeffs = None
                rhs_combo += rhs
            else:
                lhs_none, rhs = None, line.strip().lstrip()
                rhs_combo += rhs

        logging.debug(('lhs', lhs))
        logging.debug(('rhs', rhs_combo))

        reactants = lhs.split('+')
        products = rhs_combo.split('+')
        # remove trailing and leading whitespace
        reactants = [reactant.strip().lstrip().replace(' ', '') for reactant in reactants]
        products = [product.strip().lstrip().replace(' ', '') for product in products]

        logging.info(('reactants', reactants))
        logging.info(('products', products))
        logging.info(('coefficients', coeffs))

        N_reactants = len(reactants)
        if coeffs is not None:
            coeffs_list = coeffs.split(',')
            N_coeffs = len(coeffs_list)
        else:
            coeffs_list = []
            N_coeffs = 0

        equation_dict = dict()
        equation_dict['MUSICA name'] = label + ' (' + str(eqn_number) + ')'
        eqn_number += 1

        if 'hv' in lhs:
            equation_dict['type'] = 'PHOTOLYSIS'
        # MICM    A exp(C / T) (T / D)^B (1 + E p)
        # Mozart  A exp(B / T)
        elif N_coeffs == 1 and '' in products:
            equation_dict['type'] = 'FIRST_ORDER_LOSS'
            equation_dict['A'] = float(coeffs_list[0])
        elif N_coeffs == 1:
            equation_dict['type'] = 'ARRHENIUS'
            equation_dict['A'] = float(coeffs_list[0])
            equation_dict['B'] = 0
            equation_dict['C'] = 0
        elif N_coeffs == 2:
            equation_dict['type'] = 'ARRHENIUS'
            equation_dict['A'] = float(coeffs_list[0])
            equation_dict['B'] = 0
            equation_dict['C'] = float(coeffs_list[1])
        # MICM   k = A (T / 300)^B
        # Mozart k = A (300 / T)^B
        elif N_coeffs == 5:
            equation_dict['type'] = 'TROE'
            equation_dict['k0_A'] = float(coeffs_list[0])
            equation_dict['k0_B'] = - float(coeffs_list[1])
            equation_dict['kinf_A'] = float(coeffs_list[2])
            equation_dict['kinf_B'] = - float(coeffs_list[3])
            equation_dict['Fc'] = float(coeffs_list[4])
        else:
            equation_dict['type'] = 'UNKNOWN ' + equation_dict['MUSICA name']

        logging.info(equation_dict['type'])

        equation_dict['reactants'] = dict()
        equation_dict['products'] = dict()

        for reactant in reactants:
            if 'hv' in reactant:
                pass
            else:
                if len(reactant) > 0:
                    x, M = parse_term(reactant)
                    if x > 1:
                        equation_dict['reactants'][M] = {'qty': x}
                    else:
                        equation_dict['reactants'][M] = { }
        for product in products:
            if len(product) > 0:
                x, M = parse_term(product.replace('*', ' '))
                equation_dict['products'][M] = {'yield': x}

        if 'UNKNOWN' in equation_dict['type']:
            unknowns.append(equation_dict['type'])
        else:
            equations.append(equation_dict)

    return equations


if __name__ == '__main__':

    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--logfile', type=str,
        default=sys.stdout,
        help='log file (default stdout)')
    parser.add_argument('--input_dir', type=str,
        default=os.path.join('..', 'am4', 'mozart'),
        help='input config directory')
    parser.add_argument('--input_name', type=str,
        default='am4',
        help='config name')
    parser.add_argument('--micm_dir', type=str,
        default=os.path.join('..', 'am4', 'micm'),
        help='MICM output species config file')
    parser.add_argument('--mechanism', type=str,
        default='AM4',
        help='mechanism name')
    parser.add_argument('--debug', action='store_true',
        help='set logging level to debug')
    parser.add_argument('--show_json', action='store_true',
        help='show JSON output')
    args = parser.parse_args()

    """
    Setup logging
    """
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(stream=args.logfile, level=logging_level)

    """
    Read Mozart config files
    """
    lines = read_mozart_config(args.input_dir, args.input_name)

    """
    Split configs by section
    """
    sections = split_by_section(lines)
    for section in sections:
        logging.info('____ section %s ____' % section)
        for line in sections[section]:
            logging.debug(line)
        print('\n')

    """
    Generate MICM species JSON from #Fixed section
    """
    fixed_species_json = parse_species(sections['#Fixed'], fixed=True)

    """
    Generate MICM species JSON from #Variable section
    """
    variable_species_json = parse_species(sections['#Variable'])

    """
    Generate MICM equaitons JSON from #Equations section
    """
    equations_json = parse_equations(sections['#Equations'])

    """
    Assemble MICM species JSON
    """
    micm_species_json = {'camp-data': fixed_species_json + variable_species_json}
    micm_species_json_str = json.dumps(micm_species_json, indent=4)
    micm_species_yaml = yaml.dump(json.loads(micm_species_json_str))
    if args.show_json:
        logging.info('____ MICM species ____')
        logging.info(micm_species_json_str)
        print('\n')
    # logging.info(micm_species_yaml)

    """
    Assemble MICM reactions JSON
    """
    micm_reactions_json = {'camp-data':
        [{'name': args.mechanism, 'type': 'MECHANISM', 'reactions': equations_json}]}
    micm_reactions_json_str = json.dumps(micm_reactions_json, indent=4)
    micm_reactions_yaml = yaml.dump(json.loads(micm_reactions_json_str))
    if args.show_json:
       logging.info('____ MICM reactions ____')
       logging.info(micm_reactions_json_str)
       print('\n')
    # logging.info(micm_reactions_yaml)

    """
    Write MICM JSON
    """
    micm_mechanism_dir = os.path.join(args.micm_dir, args.mechanism)
    if not os.path.exists(args.micm_dir):
        os.mkdir(args.micm_dir)
    if not os.path.exists(micm_mechanism_dir):
        os.mkdir(micm_mechanism_dir)
    with open(os.path.join(micm_mechanism_dir, 'species.json'), 'w') as f:
        json.dump(micm_species_json, f, indent=4)
    with open(os.path.join(micm_mechanism_dir, 'reactions.json'), 'w') as f:
        json.dump(micm_reactions_json, f, indent=4)
    with open(os.path.join(micm_mechanism_dir, 'species.yaml'), 'w') as f:
        yaml.dump(micm_species_yaml, f, indent=4, default_flow_style=False, allow_unicode=True)
    with open(os.path.join(micm_mechanism_dir, 'reactions.yaml'), 'w') as f:
        yaml.dump(micm_reactions_yaml, f, indent=4, default_flow_style=False, allow_unicode=True)

    """
    Write unknown equation labels
    """
    with open(os.path.join(micm_mechanism_dir, 'unknowns.txt'), 'w') as f:
        for line in unknowns:
            f.write(line + "\n")

