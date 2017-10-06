from setuptools import setup, find_packages, Extension

setup(name='piperine',
    version='0.3',

    packages=['piperine', 'piperine.tests'],
    install_requires=["numpy","scipy", "stickydesign", "peppercompiler"],
    include_package_data=True,

    package_data={
         'piperine':['piperine/data/*', 'piperine/tests/test_data/*']
    },
    dependency_links=["http://www.nupack.org", "http://dna.caltech.edu/DNA_Sequence_Design_Tools/"],
    exclude_package_data={'': ['*.pyc', '*config_choi*']},
    description='Collection of software associated with the Winfree lab at Caltech',
    author='James Parkin',
    url='http://www.dna.caltech.edu/DNA_Sequence_Design_Tools/',
    license='GNU GPLv3',
    keywords='niranjan srinivas DNA dna piperine pepper winfree crn james parkin',
    zip_safe=False)
