Crispi.py is produced by the French Toast Mafia team in Professor Ian Holmes’ Bioengineering 131 class. This script intakes a fasta file(.fa) and outputs(.txt) the top ten NGG cas9 sites.

Notes:  1) Only evaluates a randomly selected SelectionSize, [NGG/NAG] sites
        2) Only looks up to MaxMismatches base pair mismatches
        3) Only outputs the top NumberOutputted sites of the selected SelectionSize into the output file
        4) This script needs the regex library to run: https://pypi.python.org/pypi/regex

Usage: python3 Crispi.py Input.fasta OutputName MaxMismatches SelectionSize NumberOutputted [NGG/NAG] (Use '-' for default)

— OutputName.txt —

-- Parameters --
Total Sequence Length: --
Max Mismatches: -- , Selection Size: --, Cas9 Site: --

Sequence (start, end) Score: ##

CrispiSeg.py is also produced by the French Toast Mafia team in Professor Ian Holmes’ Bioengineering 131 class. This script intakes a genome file(.fa), a region of that genome (.fa) and outputs(.txt) the all the specified cas9 sites.

Notes: 1) Only looks up to MaxMismatches base pair mismatches
2) Only outputs the top NumberOutputted sites of the selected SelectionSize into the output file
3) This script needs the regex library to run: https://pypi.python.org/pypi/regex

Usage: python3 CrispiSeg.py Genome.fasta Region.fasta OutputName MaxMismatches NumberOutputted [NGG/NAG]  (Use '-' for default)

— OutputName.txt —

-- Parameters --
Total Sequence Length: --
Max Mismatches: -- , Selection Size: --, Cas9 Site: --, Total Sites: --

Sequence (start, end) Score: ##