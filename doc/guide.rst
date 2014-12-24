==========
User guide
==========

Run WhatsHap like this::

	python3 -m whatshap input.vcf input.bam > phased.vcf

Phasing information is added to the VCF file in a way that is compatible with
GATK’s ReadBackedPhasing. That is, the HP tag denotes which set of phased
variants a variant belongs to.


Debugging
=========


	$ gdb python3
	(gdb) run -m nose

After you get a SIGSEGV, let gdb print a backtrace:

	(gdb) bt

