# Cytoscape analysis:

#### Step 1: get a list of unique gene families within the decision trees:
````
for i in *.dot; do grep "GeneFamily" $i | cut -f2  -d '"' | cut -f2 -d '-' > $i.txt; done
``````
#### Step 2: Get a protein sequence for each gene family

This required manually downloading the gene family protein sequences and selecting one sequence for each gene family.

#### Step 3: put proteins into String to find connections.
https://string-db.org/cgi/input?sessionId=bbEZAOFXr8HJ&input_page_active_form=COG_multiple_sequences

#### Step 4: Download network and visulise in Cytoscape

#### Step 5: Reduce the network

The network is very large, we reduced the connections by only including pairwise connections present in the network if BOTH gene families are ALSO present within at least one decision tree model.

To do this:

```
python File_compare.py | grep "both" > ampicillin_cytoscape.txt
```
where File_compare.py is:

``` 
# Open file1.txt and file2.txt
with open("String_nodes.txt") as f1, open("ampicillin_EGGNOG_NEW.arff.dot.txt") as f2:
    # Create a set of COG values from file2.txt
    file2_values = set(f2.read().splitlines())

    # Loop through each line in file1.txt
    for line in f1:
        node1, node2 = line.strip().split("\t")
        # Check if both node1 and node2 are in file2_values set
        if node1 in file2_values and node2 in file2_values:
            print(f"{node1} and {node2} are both present in file2.txt")
        else:
            print(f"{node1} and/or {node2} are not present in file2.txt")
```

This will produce a file which can be opened in microsoft excel (tab delimited). 
Paste in the *_cytoscape.txt files, the column beside the data put the antibiotic name e.g.

```
COG0038	and	COG5001	are	both	present	in	file2.txt	ciprofloxacin
COG0038	and	COG0053	are	both	present	in	file2.txt	ciprofloxacin
COG0038	and	COG3620	are	both	present	in	file2.txt	ciprofloxacin
COG0050	and	COG0154	are	both	present	in	file2.txt	ciprofloxacin
COG0110	and	COG3620	are	both	present	in	file2.txt	ciprofloxacin
COG0110	and	COG5001	are	both	present	in	file2.txt	ciprofloxacin
COG0154	and	COG5001	are	both	present	in	file2.txt	ciprofloxacin
COG0433	and	COG3620	are	both	present	in	file2.txt	ciprofloxacin
COG0433	and	COG0827	are	both	present	in	file2.txt	ciprofloxacin
COG0824	and	COG5001	are	both	present	in	file2.txt	ciprofloxacin
COG0827	and	COG3886	are	both	present	in	file2.txt	ciprofloxacin
COG1246	and	COG2916	are	both	present	in	file2.txt	ciprofloxacin
COG1246	and	COG5001	are	both	present	in	file2.txt	ciprofloxacin
```
after you have pasted each file into the file (so all the nodes and antibiotics are in the same column), delete the columns with the words 'and', 'are', 'both', 'present', 'in', and 'file2.txt'.

Then we will save this as String_intractions_by_dt.txt

#### python script to find which rows of the original network are needed:

```
# read in packages:
import pandas as pd

# read in data:
decision_tree_data = pd.read_csv("String_intractions_by_dt.txt", sep='\t')
string_interactions = pd.read_csv("string_interactions_short.tsv", sep='\t')

# join data:
combined_data = pd.merge(left=decision_tree_data, right=string_interactions, on=['#node1', 'node2'], how='left')

# output as a csv file:
combined_data.to_csv("cytoscape_input.tsv", sep='\t', index=False)
```
#### The network can be opened in cytoscape. 

In the analysis I coloured the edges via antibiotic model and the node size was proportional to the amount of models the gene family (node) was present in.




