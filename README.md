# oncoplot
1. run make_cmds.vcf2maf() with input.yaml to get a shell script to convert vcfs to mafs.
2. Use list of mafs and pass into make_cmds.merge_library_mafs() to get shell script to merge mafs into single mafs per library
3. Use list of library mafs and pass into make_cmds.annotate() to get copy number annotations for genes of interest
4. use make_cmds.merge_library_data() to merge library mafs into a single maf representing the input yaml. Use the same function to merge copynumber annotoations nto a single tabl for the input yaml.
5. use merges copy number and maf in make_oncoplot.ipynb to generate oncoplot.
