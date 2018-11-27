#!/usr/bin/python
# getBranchFileEnding.py
# Created: 2018-11-27 17:29:03.923313
# Author: Fabian Krinner
import os, sys

def getBranchFileEnding():
	inFileName = "./branchFileEnding.h"
	with open(inFileName,'r') as inFile:
		for line in inFile.readlines():
			if line.strip().startswith("//"):
				continue
			if "branchFileEnding" in line:
				return line.split()[~0].replace('"','').replace(';','')
	raise IOError("Could not load branchFileEnding from '" + inFileName + "'")

def main():
	print getBranchFileEnding()

if __name__ == "__main__":
	main()
