# Modeller loop

Pipeline for structure threading with modeller.  Requires a licensed copy of modeller (https://salilab.org/modeller/).  We use the version distributed through SBGRID (https://sbgrid.org/software/).

## Protocol

 1. Start with a multiple sequence alignment (fasta) and a structure model (pdb) that matches the region you want to perform homology threading onto.  Place both files into your working directory.  If your MSA was generated in JalView and output with trailing slashes, be sure to remove the slashes in advance using `slashRemove.py`.

 2. Run `automodel_setup.py` to prepare PIR-formatted alignment files for each of your sequences.  A new folder alled `ali_files` will be created in your working directory.

`python automodel_setup.py input_fasta_file`

 3. Adjust parameters in the `modeller_loop_new.py` file using a text editor.  You will need to specify the name of your pdb file to thread onto and an output file prefix.
 
 4. Run the modeller loop.  It should take 15-30 s per sequence you are threading.

`modeller modeller_loop_new.py`

 5. Gather the threaded files into a new subdirctory called `modeller_out`.

`python gather_modeller_outputs.py`


Issues?  Drop me a line at bell@molbio.mgh.harvard.edu.

> Written with [StackEdit](https://stackedit.io/).
