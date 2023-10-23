"""
gather_modeller_outputs.py

"""

import sys
import os
import shutil


def main():
	# Check for the required output folders, create them if they do not exist.
	neededDirs = ["modeller_out"]
	for neededDir in neededDirs:
		if os.path.exists(neededDir) == False:
			os.mkdir(neededDir)

	# Gather all files that contain ".B99990001.pdb"
	for each_file in os.listdir():
		if ".B99990001.pdb" in each_file:
			index = each_file.find(".B99990001.pdb")
			shutil.move(each_file, "modeller_out/"+each_file[:index]+".pdb")



if __name__ == '__main__':
	main()