from setuptools import setup, find_packages, Extension

from distutils.command.build import build
from setuptools.command.install import install
from setuptools.command.develop import develop

class install_with_spurious(install):
    def run(self):
        import os
        os.system("cc -Wall -O3 peppersuite/SpuriousDesign/spuriousSSM.c -o peppersuite/PepperCompiler/_spuriousSSM -lm")
        install.run(self)

class build_with_spurious(build):
    def run(self):
        import os
        os.system("cc -Wall -O3 peppersuite/SpuriousDesign/spuriousSSM.c -o peppersuite/PepperCompiler/_spuriousSSM -lm")
        build.run(self)

class develop_with_spurious(develop):
    def run(self):
        import os
        os.system("cc -Wall -O3 peppersuite/SpuriousDesign/spuriousSSM.c -o peppersuite/PepperCompiler/_spuriousSSM -lm")
        develop.run(self)

spurious_ext = Extension('peppersuite.SpuriousDesign.spuriousSSM', 
                         sources=['piperine/SpuriousDesign/spuriousSSM.c'])

setup(name='piperine',
    version='0.2',
    
    packages=['piperine',
              'piperine.SpuriousDesign', 
              'piperine.PepperCompiler',
              'piperine.PepperCompiler.design'],
    install_requires=["numpy","scipy"],
    include_package_data=True,
    # ext_modules = [spurious_ext],
    
    cmdclass={'install': install_with_spurious, 
              'build': build_with_spurious, 
              'develop': develop_with_spurious},
    
    entry_points={ 'console_scripts': [
        #'pepper-compiler = PepperCompiler.compiler:main',
        #'pepper-design-spurious = PepperCompiler.design.spurious_design:main',
        #'pepper-finish = PepperCompiler.finish:main',
        #'pepper-config = PepperCompiler.config:main',
        'spuriousSSM = piperine.PepperCompiler._spuriousSSM_wrapper:main'
        ]},
    
    scripts=['piperine/SpuriousDesign/spuriousC.py', 
             'piperine/SpuriousDesign/countspurious.py', 
             'piperine/SpuriousDesign/findspurious.py', 
             'piperine/SpuriousDesign/findspuriouswobble.py', 
             'piperine/SpuriousDesign/anneal.py',
             'piperine/PepperCompiler/compiler.py',
             'piperine/PepperCompiler/config.py',
             'piperine/PepperCompiler/finish.py',
             'piperine/PepperCompiler/design/avoid_design.py',
             'piperine/PepperCompiler/design/random_design.py',
             'piperine/PepperCompiler/design/spurious_design_constantine.py',
             'piperine/PepperCompiler/design/spurious_design.py',
             'piperine/PepperCompiler/design/new_loading.py',
             'piperine/PepperCompiler/design/PIL_parser.py',
             'piperine/PepperCompiler/quickargs.py',
             'piperine/designer.py'],
    package_data={
         'piperine':['piperine/SpuriousDesign/examples/*',
                     'piperine/SpuriousDesign/*.c',
                     'piperine/SpuriousDesign/Makefile',
                     'piperine/SpuriousDesign/make_files.zip',
                     'piperine/SpuriousDesign/make_spuriousC',
                     'piperine/data/*',
                     'piperine/PepperCompiler/doc/*',
                     'piperine/PepperCompiler/examples/*',
                     'piperine/PepperCompiler/README'],
         'piperine.SpuriousDesign':['examples/*'],
         'piperine':['data/*'],
         'piperine.PepperCompiler':['_spuriousSSM'],
         'piperine.SpuriousDesign':['*',
                                        'make_files.zip'
                                        'examples/*',
                                        'README']
    },
    dependency_links=["http://www.nupack.org", "http://dna.caltech.edu/DNA_Sequence_Design_Tools/"],
    classifiers=[
        'Programming Language :: Python :: 2.7'
    ],
    exclude_package_data={'': ['*.pyc', '*config_choi*']},
    description='Collection of software associated with the Winfree lab at Caltech',
    author='Various',
    url='http://www.dna.caltech.edu/DNA_Sequence_Design_Tools/',
    license='paypperview',
    keywords='niranjan srinivas DNA dna piperine pepper winfree crn james parkin',
    zip_safe=False)
