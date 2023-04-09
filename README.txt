This is to display our work. The code is well commented, but could benefit from a walk-through example of its use. User guide to follow if there are any requests.

-Object class can be found in "Creases_to_Input.py". 

-For examples of its use, see "Run.py"

As can be seen in "Run.py", there are two main ways to set a crease pattern. One is using a file and the other involves a command prompt dialogue. A file with example input for the water bomb can be found in "input_creases.txt". Note that equilibrium angles for the creases are being set inside the Class. This is for ease of use during our experiments.



"input_creases.txt" contains the total number of creases in the first line (7) and proceeds to define them in the next 7 blocks of lines of 4:

-First line is the mode for input (1 in all of them i.e. given end-points of crease in SCALED coordinates)
-Second line is one of the endpoints of the crease (e.g. '1/2 1/2' refers to the centre of mass of the piece of paper irrespective of its dimensions. Or '1 1' specifies top-right corner of paper)
-Third line specifies the other end-point for the crease.
-Fourth line is mountain (-1), or valley (1).

Constructed files are then to be used as prescribed by origamiMD