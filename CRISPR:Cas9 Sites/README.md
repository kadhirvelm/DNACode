Crispi.py is produced by the French Toast Mafia team in Professor Ian Holmes’ Bioengineering 131 class. This script intakes a fasta file(.fa) and outputs(.txt) the top ten NGG cas9 sites.

Notes:  1) Only evaluates a randomly selected SelectionSize, [NGG/NAG] sites
        2) Only looks up to MaxMismatches base pair mismatches
        3) Only outputs the top NumberOutputted sites of the selected SelectionSize into the output file
        4) This script needs the regex library to run: https://pypi.python.org/pypi/regex

Usage: python3 Crispi.py Input.fasta OutputName MaxMismatches SelectionSize NumberOutputted [NGG/NAG]

— OutputName.txt —
Sequence Length: ##
Sequence (start, end) Score: ##