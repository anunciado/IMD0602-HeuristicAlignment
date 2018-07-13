# Heuristic Alignment

A multiple sequence alignment is intended to perform a sequence alignment of three or more biological sequences, such as protein, DNA or RNA. In all alignments it is considered a substitution matrix, which punctuates each position of the alignment and makes the appropriate modifications (consider match / mismatch, open or extend gap). With this, homology can be inferred, which may result in a phylogenetic analysis to evaluate the evolutionary origins shared by the sequences.

However, for the alignment of multiple sequences of biologically relevant length one can take a significant time, so it becomes necessary computational algorithms that can optimize the execution time to analyze the alignments. Most multi-sequence alignment programs use heuristic methods, which try to identify good alignment in all runs, often failing to find the best alignment possible.

In this context, the project receives a 'query.fa' file containing the references sequences and a 'database.fa' file so that they contain a database of sequences to be aligned with the reference sequences. First it will be called the readQuery method that will put its contents in two lists, one with the name of the sequences and another with the sequences, with these lists will be called the method readDatabase that will iterate over the lists passed with parameter to the method, where each interaction will be placed in a hash table as a key of a size 11 word drawn from the reference sequence and as value its position in the original sequence. Soon after, each word in the sequence of the database will be searched and stored in a dictionary, the positions of the sequences found, to recognize the two largest set of sequence matches of 11 bases of the reference sequence contained in the database sequence. Then, each iteration will align the reference sequences with a database sequence through the localAlignment function, which uses the score function that gives two base pairs, returns an alignment value, which will do the local alignment between the positions of the two largest match sequences and will return an alignment and score for that alignment.

At the end of the loop, the identity method will be called that will store the values in a dictionary, where the keys are the GenBank codes of the reference sequences and the values of each key will be lists containing useful information (number of matchs, sequence name database, sequence reference name, score, sequence alignment, reference sequence size, database sequence size), and in the end, returns the 10 best alignments based on the score value of the alignments. Subsequently, the final dictionary with the 10 best alignments for each reference sequence in 'query.fa' will be shown in the terminal.

### Prerequisites

You will need to install the modules below to run the program: 
* [python 3.7 or greater](https://www.python.org/downloads/release/python-370/)

### Running

There are two ways to run the program:

* Compile the IDE (PyCharm - Python IDE):
1. Just open the IDE
2. Import the project folder as a Project
3. Choose Run scheduling on the context menu.
4. From this it only interacts with the system and add in script parameters box contents:

* Compile by terminal:
1. Enter the src folder and run the following command:
```
python heuristicAlignment.py
```
2. From this it only interacts with the system.

## Built With

* [PyCharm](https://www.jetbrains.com/pycharm/) - A IDE used

## Authors
### Developers: 
* **Luís Eduardo Anunciado Silva ([cruxiu@ufrn.edu.br](mailto:cruxiu@ufrn.edu.br))** 
### Project Advisors: 
* **Cesar Renno Costa ([cesar@imd.ufrn.br](mailto:cesar@imd.ufrn.br))**
* **Jorge Estefano Santana De Souza ([jorge@imd.ufrn.br](mailto:jorge@imd.ufrn.br))**

See also the list of [contributors](https://github.com/cruxiu/IMD0602-HeuristicAlignment/contributors) who participated in this project.

## License

This project is licensed under the GPL 3.0 - see the [LICENSE](LICENSE) file for details
