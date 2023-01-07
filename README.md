# Context
Janina has a simple DNA design problem that is suitable for automation. I am attempting to write that automation. This directory stores data and code towards those ends.

# Specs
Given 
1. excel spreadsheet of DNA sequences representing a CDS (ATG -> STOP)
2. a species specific degron sequence (protein space) to replace

Do
	1. Append a static prefix and suffix
	2. Search and replace all instances of degron with
		1. deletion (cut nucleotide sequence) 
		2. poly Alanine (replace nucleotides with same alanine codon)
	3. return two excel spreadsheets with the replaced sequences

# Caveats
degron sequence will change between runs! 

some sequences in the input spreadsheet will not contain any instances of the degron
