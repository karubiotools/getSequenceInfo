# ![Logo getSequenceInfo](/logo_getSequenceInfo.png)

A simple [Perl](https://en.wikipedia.org/wiki/Perl) script allowing to get sequence information from GenBank, RefSeq or ENA sequence repositories.

# Requirements
Perl (version 5.26 or greater) must be available in your system to run getSequenceInfo. If your Operating System (OS) is Windows, you can get Perl by installing [Strawberry Perl](http://strawberryperl.com/). If necessary, please see information on [how to launch](https://www.digitalcitizen.life/7-ways-launch-command-prompt-windows-7-windows-8) or [how to use](https://www.digitalcitizen.life/command-prompt-how-use-basic-commands) the Command Prompt in Windows.
When using Unix OS (Linux or Mac), Perl is generally already installed. But if it is not the case, you can see [this page](https://learn.perl.org/installing/unix_linux.html) for its installation. You can follow this [wiki page](https://en.wikibooks.org/wiki/Guide_to_Unix/Explanations/Shell_Prompt) for information about the Shell Prompt.
You can then check the installation by typing the following command:
```bash
perl -v
```

# Quick Installation
Please first verify that Perl is installed in your system by following the above requirments.
You probably need to install the X11 development package first.
On Debian or Ubuntu, this is the package libx11-dev: ```sudo apt-get install libx11-dev```  
On CentOS, RedHat, or Fedora, this is the package libX11-devel.  
MacOS users may need Xcode/XQuartz and Fink programs.
## Linux or MacOS (Unix)
```bash
git clone https://github.com/dcouvin/getSequenceInfo.git
cd getSequenceInfo/
bash install/installer_Unix.sh
```

## Windows
Users can also install the tool by running the installer_Windows.bat file (double-click)
```bash
install\installer_Windows.bat
```

# How to use
The tool can be used directly with the command line or using a graphical user interface (GUI).
## GUI version
The user can launch the GUI version of the tool (getSequenceInfoGUI.pl) either by executing it (double click) or by typing the following command:
```bash
perl getSequenceInfoGUI.pl
```
## Command line version
We can type the following command to display the help message:
```bash
perl getSequenceInfo.pl -h
```
Help message:
```
	Name: 
		getSequenceInfo.pl
	
	Synopsis:
		A Perl script allowing to get sequence information from GenBank RefSeq or ENA repositories.
		
	Usage:
	  perl getSequenceInfo.pl [options]
	  examples: 
	     perl getSequenceInfo.pl -k bacteria -s "Helicobacter pylori" -le "Complete Genome" -date 2019-06-01 
	     perl getSequenceInfo.pl -k viruses -n 5 -date 2019-06-01
	     perl getSequenceInfo.pl -k "bacteria" -taxid 9,24 -n 10 -c plasmid -dir genbank -o Results
	     perl getSequenceInfo.pl -ena BN000065
	     perl getSequenceInfo.pl -fastq ERR818002
	     perl getSequenceInfo.pl -fastq ERR818002,ERR818004
						 	
	Kingdoms:
		archaea
		bacteria
		fungi
		invertebrate
		plant
		protozoa
		vertebrate_mammalian
		vertebrate_other
		viral
	
	Assembly levels:
		"Complete Genome"
		Chromosome
		Scaffold
		Contig 
	
	General:
		-help or -h			displays this help 	
		-version or -v			displays the current version of the program
		
	Options ([XXX] represents the expected value):
		-directory or -dir [XXX]	allows to indicate the NCBIs nucleotide sequences repository (default: genbank)
		-get or -getSummaries [XXX]	allows to obtain a new assembly summary files in function of given kingdoms (bacteria,fungi,protozoa...)	
		-k or -kingdom [XXX]		allows to indicate kingdom of the organism (see the examples above)
		-s or -species [XXX]		allows to indicate the species (must be combined with -k option)
		-taxid [XXX]			allows to indicate a specific taxid (must be combined with -k option)
		-assembly_or_project [XXX]	allows to indicate a specific assembly accession or bioproject (must be combined with -k option)
		-date [XXX]			indicates the release date (with format yyyy-mm-dd) from which sequence information are available
		-le or -level [XXX]		allows to select a specific assembly level (e.g. "Complete Genome")
		-o or -output [XXX]		allows users to name the output result folder
		-n or -number [XXX]		allows to limit the total number of assemblies to be downloaded
		-c or -components [XXX]		allows to select specific components of the assembly (e.g. plasmid, chromosome, ...)
		-ena [XXX] 			allows to download report and fasta file given a ENA sequence ID 
		-fastq [XXX]			allows to download FASTQ sequences from ENA given a run accession (https://ena-docs.readthedocs.io/en/latest/faq/archive-generated-files.html)
		-log				allows to create a log file
```
## Singularity container
This tool can also be used from a [singularity](https://sylabs.io/singularity/) container. Make sure singularity is installed on your machine, then run the following commands:
```bash
sudo singularity build getSequenceInfo.simg getSequenceInfo.def
```
Example of use:
```bash
singularity exec -B $PWD getSequenceInfo.simg perl /usr/local/getSequenceInfo.pl -h
```
