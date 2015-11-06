The PySeqSimulator.py file does the following:

Use '--calc fasta.fasta' to create a random sequence based on fasta.fasta
Or Use '--load params' to create a random sequence based on the parameters in params
Follow either command above with '--save file' to save the logistics data to file.

The ReverseCompliment file does the following:

Run it with a fasta file and it'll write the reverse compliment to the terminal.

The InformationContent file does the following:

Run it with a Stockholm formatted file as its input and it'll spit out the top ten 
lowest Shannon entropy columns, and the top fifty mutual information columns.