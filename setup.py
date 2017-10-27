from setuptools import setup, find_packages, Extension

setup(name='piperine',
    version='0.4a',

    packages=['piperine', 'piperine.tests', 'piperine.Srinivas2017'],
    install_requires=["numpy","scipy", "stickydesign", "peppercompiler"],
    include_package_data=True,

    package_data={
         'piperine':['piperine/data/*', 'piperine/tests/test_data/*']
    },
    dependency_links=["http://www.nupack.org", "http://dna.caltech.edu/DNA_Sequence_Design_Tools/"],
    entry_points={ 'console_scripts': [
        'piperine-design = piperine.designer:main',
        'piperine-score = piperine.designer:score']},
    exclude_package_data={'': ['*.pyc', '*config_choi*']},
    description='A pipeline to generate DNA sequences that implement abstract CRNs.',
    author='James Parkin',
    url='http://www.dna.caltech.edu/DNA_Sequence_Design_Tools/',
    license='GNU GPLv3',
    keywords='niranjan srinivas DNA dna piperine pepper winfree crn james parkin',
    zip_safe=False)
