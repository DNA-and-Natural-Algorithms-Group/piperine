from setuptools import setup, find_packages, Extension

setup(name='piperine',
    version='0.4a',

    packages=['piperine', 'piperine.tests', 'piperine.Srinivas2017', 'piperine.Chen2013'],
    install_requires=["numpy","scipy","stickydesign >=0.8.4,!=0.9.0.a1","peppercompiler >= 0.1.4"],
    include_package_data=True,

    package_data={
         'piperine':['data/*', 'tests/test_data/*']
    },
    dependency_links=["http://www.nupack.org", "http://dna.caltech.edu/DNA_Sequence_Design_Tools/"],
    entry_points={ 'console_scripts': [
        'piperine-design = piperine.commandline_utilities:design',
        'piperine-select = piperine.commandline_utilities:select',
        'piperine-score = piperine.commandline_utilities:score']},
    exclude_package_data={'': ['*.pyc', '*config_choi*']},
    description='A pipeline to generate DNA sequences that implement abstract CRNs.',
    author='James Parkin',
    url='http://www.dna.caltech.edu/DNA_Sequence_Design_Tools/',
    license='GNU GPLv3',
    keywords='niranjan srinivas DNA dna piperine pepper erik winfree crn james parkin',
    zip_safe=False)
