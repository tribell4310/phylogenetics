from modeller import *
from modeller.automodel import *
import os


# USER adjusts these parameters
my_pdb_file = "my_pdb.pdb"
my_prefix_name = "isocitrate-synthase"


# Loop through each input for human structure template
for each_file in os.listdir("./ali_files"):
	# Generate alignment files
	env = Environ()
	aln = Alignment(env)
	mdl = Model(env, file=my_pdb_file, model_segment=('FIRST:A','LAST:A'))
	aln.append_model(mdl, align_codes='syntheticSeq', atom_files=my_pdb_file)
	aln.append(file="./ali_files/"+each_file)
	aln.align2d(max_gap_length=50)
	aln.write(file=my_prefix_name+'-'+each_file, alignment_format='PIR')

	# Structure threading
	a = AutoModel(env, alnfile=my_prefix_name+'-'+each_file, knowns="syntheticSeq", sequence=each_file[:-4], assess_methods=(assess.DOPE, assess.GA341))
	a.starting_model = 1
	a.ending_model = 1
	a.make()

