# Piperine
Piperine is a software tool for the automated design of DNA molecules that mimic abstract chemical reaction networks (CRNs) in a test tube.
It combines sequence design tools created by the DNA and Natural Algorithms Group with heuristic sequence quality measures.
These heuristic functions detect sequence motifs that are detrimental to DNA strand displacement reactions; avoiding these motifs promotes fidelity of the _in vitro_ DNA reaction to the target CRN.

Major Contributors :

* James Parkin

* Niranjan Srinivas

* Erik Winfree

* Harel Dor

Thanks to : 

* Chris Thachuk 

* Constantine Evans 

Basic use of Piperine involves providing CRNs specified in text files and receiving a set of candidate DNA implementations.
From this list of implementations, Piperine will select one as the optimal candidate relative to the set.
Users can use command line utilities to adjust the design constraints and number of candidate implementations generated.
Advance use may involve writing custom packages that define additional approaches to emulating CRNs with DNA strand displacement cascades.

## Installation

### OS requirements
Piperine has only been tested on Ubuntu and MacOSX Sierra.
Windows is unsupported for now.

### Python version
Piperine works with Python >= 2.7

### Software dependencies
Piperine depends on

* Numpy

* Scipy

* NUPACK (version 3.0.x)

* stickydesign

* peppercompiler

Make sure these are installed before installing Piperine.

### Installing dependencies
Numpy and Scipy can both be installed easily through pip or conda.
NUPACK can be [downloaded here](http://www.nupack.org/).
Stickydesign and peppercompiler can both be found [on the DNA and Natural Algorithms Group Github page](https://github.com/DNA-and-Natural-Algorithms-Group).
Be sure to install NUPACK version 3.0.x, Piperine is only compatible with this legacy version for now.
Once NUPACK is installed, make a terminal variable NUPACKHOME that points to installation destination (might be `~/Downloads/nupack3.0.6`). Do this with the command `export NUPACKHOME=`_the path to nupack_.

### Installing Piperine
For now, you can only install Piperine from the source files.
Clone or download this repository, start a terminal session, and navigate to the repository folder.
From there, use __pip__ to install Piperine using the in-place flag __-e__. 

```
>user@laptop:~/git$ git clone git@github.com:DNA-and-Natural-Algorithms-Group/piperine.git
>Cloning into 'piperine'...
>user@laptop:~/git$ ls piperine # check the contents of the piperine directory
>LICENSE piperine README.md setup.py
>
>user@laptop:~/git$ pip install -e piperine
```

To install updates, return to the git repository and execute

`git pull`

## Brief tutorial
Piperine reads in plaintext files that specify an abstract CRN, converts the CRN into a domain-level DNA specification, then designs sequences implementing the domain-level constraints.
Multiple sets of sequences are produced during each execution.

### Command line utilities
#### Sequence design and selection
The command line utility `piperine-design` requires one argument, a file specifying an abstract CRN.
An example CRN might be the following, saved in a file called `my.crn`

```
 -> B
B -> A
A + A -> A
```

This CRN can then be compiled into DNA sequences using the following call in a terminal session:

`piperine-design my.crn -n 3`

The `-n 3` argument tells Piperine to generate three candidate sets of sequences.
This execution will generate a number of intermediate files suffixed by .sys, .pil, .fixed, and .mfe.
Most users will only be interested in the candidate comparison report suffixed by score_report.txt and the files suffixed by .seqs which contain the sequence identities of the strands.
At the end of the process, the winning candidate will be announced.
Sequence sets are indexed starting at 0 and saved to filenames `my`__i__`.crn`, for index __i__.

Users may also adjust toehold generation parameters through command line arguments.
These optional arguments define the target binding energy, allowable deviation from the target energy, and maximum relative energy of unintended toehold interactions.
To override the default values, apply arguments as in :

`piperine-design my.crn --energy 7.5 --deviation 0.5 --maxspurious 0.5`

This may be necessary when generating sequences for CRNs including more than 10 toeholds.
If the default or user-specified toehold energetics parameters cannot be satisfied by the sequence designer, Piperine will suggest values that will generate a sufficient number of toeholds.

More complete information on the `piperine-design` utility may be found in its documentation:

`piperine-design --help`

#### Translation schemes
A translation scheme is a set of rules that define how a general CRN may be translated into a DNA implementation.
For Piperine, translation schemes are Python packages that provide the information and functions required for Piperine to design and score DNA sequences.
Two translation schemes are included with Piperine: Chen2013 and Srinivas2017.
For more information on these schemes, review the publications listed in the references below.

By default Piperine will use the Srinivas2017 scheme when designing sequences.
To use the Chen2013 scheme, use the `--translation_scheme` argument.

`piperine-design my.crn --translation_scheme Chen2013`

## References
1. Chen, Y. J. et al. Programmable chemical controllers made from DNA. Nat. Nanotechnol. 8, 755–762 (2013).
