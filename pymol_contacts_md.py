#! /usr/bin/env python3

"""
Created on 9 Oct. 2023
"""

import argparse
import logging
import os
import re
import sys

import pandas
import pymol

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.1.0"


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[args.log_level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def extract_roi(roi_to_extract):
    """
    Extract the region of interest (roi) start's and stop's coordinates.

    :param roi_to_extract: the coordinates ot the region of interest, as 100-200 i.e.,
    :type roi_to_extract: str
    :raises ArgumentTypeError: is not between 0.0 and 100.0
    :return: the region of interest (roi) start's and stop's coordinates.
    :rtype: list
    """
    pattern_roi_to_extract = re.compile("(\\d+)-(\\d+)")
    match_roi_to_extract = pattern_roi_to_extract.search(roi_to_extract)
    if match_roi_to_extract:
        pos1 = int(match_roi_to_extract.group(1))
        pos2 = int(match_roi_to_extract.group(2))
        if pos1 <= pos2:
            roi_extracted = [int(match_roi_to_extract.group(1)), int(match_roi_to_extract.group(2))]
        else:
            raise argparse.ArgumentTypeError(f"'{roi_to_extract}' argument for the option --roi is malformed, the "
                                             f"first position ({pos1}) is > to the second position ({pos2}).")
    else:
        raise argparse.ArgumentTypeError(f"'{roi_to_extract}' argument for the option --roi is malformed, it should be "
                                         f"two integers separated by an hyphen, i.e: '100-200'.")
    return roi_extracted


def create_contact(atoms_contacts):
    """
    Create the contact between the two atoms of the residues.

    :param atoms_contacts: The atoms' contacts.
    :type atoms_contacts: str
    :return: the number of contacts added.
    :rtype: int
    """
    pattern = re.compile("\\D{3}(\\d+)_(\\S+?)-\\D{3}(\\d+)_(.+)")
    # if the atom 2 belongs to a hydrogen bond, the notation will be atom2-hydrogen_donor
    pattern_atom2 = re.compile("([A-Z0-9]+)-([A-Z0-9]+)")
    number_contacts = 0
    for atom_contact in atoms_contacts.split(" | "):
        match = pattern.search(atom_contact)
        if match:
            resi_1 = match.group(1)
            atom_1 = match.group(2)
            resi_2 = match.group(3)
            atom_2_tmp = match.group(4)
            hydrogen_atom_from_hbond = None
            match_atom2 = pattern_atom2.search(atom_2_tmp)
            if match_atom2:
                # atom 2 from hydrogen bonds analysis
                atom_2 = match_atom2.group(1)
                hydrogen_atom_from_hbond = match_atom2.group(2)
            else:
                # atom 2 not from hydrogen bonds analysis
                atom_2 = atom_2_tmp
        else:
            raise(Exception(f"No match in {atom_contact} for '{pattern.pattern}'."))
        downgraded_contact = False
        p1 = pymol.cmd.select("p1", f"(resi {resi_1} and name {atom_1})")
        p2 = pymol.cmd.select("p2", f"(resi {resi_2} and name {atom_2})")
        if p1 != 1:
            raise Exception(f"create_contact function on {atom_contact}: PyMol command for p1 selection on "
                            f"{resi_1}_{atom_1} failed (select resi {resi_1} and name {atom_1}).")
        if p2 != 1:
            if hydrogen_atom_from_hbond:
                downgraded_contact = True
                logging.warning(f"create_contact function on {atom_contact}: PyMol command for p2 selection on "
                                f"{resi_2}_{atom_2} failed (select resi {resi_2} and name {atom_2}), trying on "
                                f"{resi_2}_{hydrogen_atom_from_hbond}.")
                p2_sel_str = f"(resi {resi_2} and name {hydrogen_atom_from_hbond})"
                p2 = pymol.cmd.select("p2", p2_sel_str)
                if p2 != 1:
                    raise Exception(f"create_contact function on {atom_contact}: PyMol commands for p2 selection "
                                    f"failed:\nselect resi {resi_2} and name {hydrogen_atom_from_hbond}\nselect resi "
                                    f"{resi_2} and name {hydrogen_atom_from_hbond}")
            else:
                raise Exception(f"create_contact function on {atom_contact}: PyMol commands for p2 selection "
                                f"failed:\nselect resi {resi_2} and name {atom_2}.")
        pymol.cmd.distance(atom_contact, "p1", "p2")
        pymol.cmd.hide("labels", atom_contact)
        if downgraded_contact:
            pymol.cmd.color("red", atom_contact)
            logging.warning(f"{atom_contact}: Downgraded contact between \"{resi_1}_{atom_1} and "
                            f"{resi_2}_{hydrogen_atom_from_hbond}\" instead of \"{resi_1}_{atom_1} and "
                            f"{resi_2}_{atom_2}\" set to red color instead of yellow.")
        pymol.cmd.delete("p1")
        pymol.cmd.delete("p2")
        number_contacts += 1
    return number_contacts


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or
    implied.

    Add to a structure file of a protein the contacts discovered by the plot_hbonds.py 
    (https://github.com/njeanne/plot_hbonds) or by the plot_neighbors.py scripts 
    (https://github.com/njeanne/plot_neighbors).
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--prefix", required=True, type=str,
                        help="the prefix of the path to the PyMol file. A \".pse\" extension will be added to this "
                             "prefix.")
    parser.add_argument("-s", "--structure", required=True, type=str,
                        help="the protein structure file, it can be \".pse\", '.pdb'.")
    parser.add_argument("-r", "--roi", required=False, type=str,
                        help="the Region Of Interest coordinates, the  \"residue 1 position\" must belong to this "
                             "interval to be added as a contact. The format must be two digits separated by an "
                             "hyphen, i.e: '100-200'.")
    parser.add_argument("-d", "--domains", required=False, type=str, default="",
                        help="the path to the CSV file to annotate the protein domains. The domains file is a comma "
                             "separated file, with a column named 'domain' which contains the domains names, two  "
                             "columns respectively called 'start' and 'end' with the 1-indexed coordinates of each "
                             "domain. And finally a column called 'pymol color' with the name of the color as "
                             "described here: https://pymolwiki.org/index.php/Color_Values")
    parser.add_argument("-e", "--exclude-domains", required=False, nargs="+",
                        help="the list of domains to exclude in the contacts. A list of domains names as they appear "
                             "in the \"residue 2 domain\". The arguments must be separated by spaces and if their "
                             "names contains spaces, they must be surrounded by '\"', in example: PPR \"X domain\".")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help=("the path for the log file. If this option is  skipped, the log file is created in the "
                              "output directory."))
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If this option is skipped, the log level is INFO.")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="path to the contacts CSV file from the plot_contacts.py script.")
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.prefix))
    os.makedirs(out_dir, exist_ok=True)

    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(out_dir, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")

    # get the Region Of Interest if specified
    roi = None
    if args.roi:
        try:
            roi = extract_roi(args.roi)
        except argparse.ArgumentTypeError as exc:
            logging.error(exc)
            sys.exit(1)

    contacts = pandas.read_csv(args.input, sep=",")
    pymol.cmd.load(args.structure)
    # set the protein to white color
    pymol.cmd.color("white", "all")

    domains_data = None
    if args.domains:
        logging.info(f"Coloring the protein domains using the file: {args.domains}")
        try:
            domains_data = pandas.read_csv(args.domains)
            for _, row in domains_data.iterrows():
                pymol.cmd.select(row["domain"].replace(" ", "_"), f"resi {row['start']}:{row['end']}")
                pymol.cmd.color(row["pymol color"], f"resi {row['start']}:{row['end']}")
        except FileNotFoundError as exc:
            logging.error(exc)
            sys.exit(1)

    excluded_domains_in_contacts = None
    domains = [x.lower() for x in domains_data["domain"].tolist()]
    idx_domains_error = []
    if args.exclude_domains:
        excluded_domains_in_contacts = [item.strip().lower() for item in args.exclude_domains]
        for idx in range(len(excluded_domains_in_contacts)):
            if excluded_domains_in_contacts[idx] not in domains:
                idx_domains_error.append(idx)
    for idx in reversed(idx_domains_error):
        del excluded_domains_in_contacts[idx]
        logging.warning(f"Ignored excluded domain \"{args.exclude_domains.pop(idx)}\". The domain is not present in "
                        f"the file --domains {args.domains}")

    nb_initial_contacts = 0
    nb_validated_contacts = 0
    nb_validated_residues_pairs = 0
    nb_excluded_contacts = 0
    nb_out_roi = 0
    for _, row in contacts.iterrows():
        nb_initial_contacts += row["number atoms contacts"]
        if roi and (int(row["residue 1 position"]) < roi[0] or int(row["residue 1 position"]) > roi[1]):
            nb_out_roi += row["number atoms contacts"]
            logging.debug(f"{row['residues contact']}: {row['number atoms contacts']} contacts excluded because the "
                          f"residue 1 position ({row['residue 1 position']}) is outside the Region Of Interest "
                          f"limits: {args.roi}.")
            continue
        if excluded_domains_in_contacts and row["residue 2 domain"].lower() in excluded_domains_in_contacts:
            nb_excluded_contacts += row["number atoms contacts"]
            logging.debug(f"{row['residues contact']}: {row['number atoms contacts']} contacts excluded because the second "
                          f"partner domain ({row['residue 2 domain']}) is in the list of the excluded domains.")
            continue
        # change the representation of the two residues which atoms are in contact to licorice
        logging.debug(f"contact added: {row}")
        pymol.cmd.select("tmp", f"resi {row['residue 1 position']} or resi {row['residue 2 position']}")
        pymol.cmd.show(representation="licorice", selection="tmp")
        pymol.cmd.delete("tmp")
        # create the contact
        try:
            nb_validated_contacts += create_contact(row["atoms contacts"])
        except Exception as exc:
            logging.error(exc, exc_info=True)
            sys.exit(1)
        nb_validated_residues_pairs += 1
    if nb_out_roi != 0:
        logging.warning(f"{nb_out_roi} contacts excluded because the residue 1 position was outside of the Region Of "
                        f"Interest limits: {args.roi}")
    if nb_excluded_contacts != 0:
        logging.warning(f"{nb_excluded_contacts} contacts excluded because the residue 2 domain was one of the "
                        f"excluded domains: {', '.join(excluded_domains_in_contacts)}.")
    logging.info(f"{nb_validated_contacts}/{nb_initial_contacts} contacts added for "
                 f"{nb_validated_residues_pairs}/{len(contacts)} pairs of residues.")
    out = f"{os.path.abspath(args.prefix)}.pse"
    pymol.cmd.save(out)
    logging.info(f"PyMol file written: {out}")
