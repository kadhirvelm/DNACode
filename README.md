The PySeqSimulator.py file does the following:

Use '--calc fasta.fasta' to create a random sequence based on fasta.fasta
Or Use '--load params' to create a random sequence based on the parameters in params
Follow either command above with '--save file' to save the logistics data to file.

The ReverseCompliment file does the following:

Run it with a fasta file and it'll write the reverse compliment to the terminal.

The InformationContent file does the following:

Run it with a Stockholm formatted file as its input and it'll spit out the top ten 
lowest Shannon entropy columns, and the top fifty mutual information columns.

The Pairwise Alignment file does the following:

Run it with a similarity matrix as the first argument and a sequence file as the second
to get the first two sequences inside the sequence file aligned and printed to the 
terminal.

Crispi.py file does the following:

Run it with any fasta file and it'll create an output.txt with the best CRISPR/Cas9 NGG sites
present in the file.'