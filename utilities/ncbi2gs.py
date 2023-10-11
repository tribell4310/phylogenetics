"""

10 September 2020
py3

Rename NCBI nomenclature to G_s nomenclature.

"""

import sys


def main(inFasta):
	f = open(inFasta, "r")
	lines = f.readlines()
	
	collateLines = []
	suspiciousLines = []
	for i in range(0, len(lines)):
		if lines[i][0] == ">":
			forBrackLoc = lines[i].find("[")
			revBrackLoc = lines[i].find("]")
			if forBrackLoc != -1:
				gSpaceS = lines[i][forBrackLoc+1:revBrackLoc]
				spaceLoc = gSpaceS.find(" ")
				secondSpaceLoc = gSpaceS.find(" ", spaceLoc+1)
				if secondSpaceLoc == -1:
					finalStr = gSpaceS[0:spaceLoc]+"_"+gSpaceS[spaceLoc+1:]
				else:
					finalStr = gSpaceS[0:spaceLoc]+"_"+gSpaceS[spaceLoc+1:secondSpaceLoc]
				collateLines.append(">"+finalStr+"\n")
				if finalStr[0].isupper() == False:
					suspiciousLines.append(finalStr)
			else:
				collateLines.append(lines[i])
		else:
			collateLines.append(lines[i])

	# Alert user
	if len(suspiciousLines) != 0:
		print("\nThe following suspicious lines were found, please take note:\n")
		for item in suspiciousLines:
			print(item)

	# Write out
	g = open(inFasta[:-6]+"_ncbi2gs.fasta", "w")
	for line in collateLines:
		g.write(line)
	g.write("\n")

	# Clean up
	g.close()
	f.close()


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()