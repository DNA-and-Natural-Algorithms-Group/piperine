import sys
import os
import pkg_resources

def main():
    binary_path = pkg_resources.resource_filename(__name__,'_spuriousSSM')

    os.system(binary_path + " " + " ".join(sys.argv[1:]))
