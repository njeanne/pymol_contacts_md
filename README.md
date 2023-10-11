# Purpose

From the contacts retrieved from a molecular dynamics simulation with the [plot_contacts.py script](https://github.com/njeanne/plot_contacts), a
representation of this contacts will be performed on a structure file of a protein. The `.pse` file can be visualized
with [PyMol](https://github.com/schrodinger/pymol-open-source).

# Environment

A [Conda](https://docs.conda.io/en/latest/) environment is provided in the `env` directory and can be loaded with the
command:
```shell
conda conda_env create -f conda_env/pymol_contacts_md_md_env.yml
```

# Usage

The example files provided in the `data` directory can be used to test the script.

## Basic usage

The basic usage is:
```shell
conda activate pymol_contacts_md

./pymol_contacts_md.py --prefix results/contacts --structure data/JQ679013_RPS17_ORF1_0.pse \
data/outliers_JQ679013_RPS17_ORF1.csv

conda deactivate
```
This will produce a the `results/contacts.pse` file that can be opened with [PyMol](https://github.com/schrodinger/pymol-open-source),
adding 81 contacts from 55 pairs of residues.

The persistent contacts during the Molecular Dynamics simulation between the atoms of different residues are added:

![contacts on the structure file](doc/_static/basic.png)

## Excluding domains

In the CSV input file, the domains of the second atom can be excluded, column `second partner domain`. If we want to
exclude the Hinge and the X domain, the command will be:
```shell
conda activate pymol_contacts_md

./pymol_contacts_md.py --prefix results/contacts_no-hinge-xdomain --exclude-domains Hinge "X domain" \
--structure data/JQ679013_RPS17_ORF1_0.pse data/outliers_JQ679013_RPS17_ORF1.csv

conda deactivate
```

This will produce a the `results/contacts_no-hinge-xdomain.pse` file that can be opened with
[PyMol](https://github.com/schrodinger/pymol-open-source), removing 40 contacts being parts of the excluded domains, to
finally add 41/81 contacts added for 31/55 pairs of residues.

## Selecting a Region Of Interest (ROI)

A region of interest in the CSV input file can be selected to add only the contacts of the first partner of the contact
which coordinates are in the targeted region of interest:

Here, we target an insertion which coordinates are amino acid 748 to 806:
```shell
conda activate pymol_contacts_md

./pymol_contacts_md.py --prefix results/contacts_roi-748-806 --roi 748-806 \
--structure data/JQ679013_RPS17_ORF1_0.pse data/outliers_JQ679013_RPS17_ORF1.csv

conda deactivate
```

This will produce a the `results/contacts_roi-748-806.pse` file that can be opened with
[PyMol](https://github.com/schrodinger/pymol-open-source), removing 45 contacts which first partner atom is outside the
region of interest 748 to 806, to finally adding 36/81 contacts added for 22/55 pairs of residues:

![contacts on the structure file](doc/_static/roi-748-806.png)
