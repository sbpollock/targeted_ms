# Mass_Spec_Targeting

This is a Shiny dashboard designed to determine m/z search values for one or more peptide targets, then extract retention time and intensity data against MS1 spectra using .raw files as an input. You can find the dockerhub page at https://hub.docker.com/r/sbpollock/targeted_ms.

You can use the dashboard to, for instance, monitor retention times and intensity values of control BSA peptides across runs for quality control.

Alternatively, if the retention time of your target peptides of interest is known, you can search for their presence in an experimental run with a bit better sensitivity (and quantitativeness, in my experience) than using a database search program like PEAKS.

If you need sample .raw files to test it out, https://www.ebi.ac.uk/pride/archive/projects/PXD009542 has really nice, well annotated files.