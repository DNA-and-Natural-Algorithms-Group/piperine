from setuptools import setup, find_packages
setup(
    name="piperine",
    version="0.1",
    packages=find_packages(),#exclude=['tests']),
    scripts=['piperine/designer.py'],

    install_requires=['numpy'],
    dependency_links=["http://www.nupack.org", "http://dna.caltech.edu/DNA_Sequence_Design_Tools/"],

    include_package_data = True,
    package_data={
        '': ['*.csv', '*.comp', '*.crn', 'tests/test_data/*'],
        'piperine': ['data/*.csv', 'data/*.comp', 'data/*.crn'],
    },
    
   # data_files = [
   #               ('comps', ['data/bimrxn.comp', 'data/leakless_and.comp', 'data/leakless_translate.comp', 
   #                          'data/leakless_flux.comp']),
   #               ('crns', ['data/oscillator.crn', 'data/Rossler.crn', 'data/RW.crn', 'data/small.crn']),
   #               ('params', ['data/dnadangle35.csv', 'data/dnadangle.csv', 'data/dnastackingbig.csv']),
   #                 ],

    # metadata for upload to PyPI
    author="James Parkin",
    author_email="jparkin@caltech.edu",
    description="Automatic design and scoring of DNA sequences for dynamical strand displacement circuits",
    license="MIT",
    keywords="dna ssa dsd crn",
    url="http://www.dna.caltech.edu/DNA_Sequence_Design_Tools/",
)
