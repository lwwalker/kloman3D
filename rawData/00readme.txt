rawData:
(1/22/13)
Contains original datafiles, sent by Sergei to Lorne 9/30/10

"Hi Lorne,

I am attaching some data on 3 monoclonal nAbs.

2F5, 447-52D and B12

There are three files for each mAB.

1). .txt is the file that has an accession number, sequence, mAB ID and the titer (lower is more sensitive)
2). .fas file are just the sequences
3). .py files are fractional coordinate files mapped by HyPhy (you can call execfile on them from your Python script)

These files have a dict (keyed on the accession number).
Each entry of the dict is a dict itself, now keyed on the position in the corresponding sequence.
Each position is a dict with two attributes:

'AA' - the amino-acid at that position and
'C' - the coordinate of that position in HXB2 system.

Sergei"

(1/23/13)
testFas.py - quick script for checking for duplicate sequences
Note: b12.fas contains duplicate sequences for 21 ID's, but all sequences for a given ID are identical
Note: 2F5.fas contains duplicate sequences for 21 ID's, but all sequences for a given ID are identical
447-52D contains no duplicate sequencessu

(1/28/13)
Added dta for 2G12, 4E10, PG9, PG16
