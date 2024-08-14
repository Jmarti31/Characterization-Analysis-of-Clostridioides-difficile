""" addKEGGPathways.py script: Identifies pathway maps associated to predicted proteins."""

import argparse
import requests
import re


def main():
    """Main Business Logic"""
    # Get arguements from the user
    args = get_args()
    input_file = args.input_file
    threshold = args.threshold
    output_file = args.output_file

    # Loads up the availble pathways in KEGG, and stroes for late usage.
    kegg_paths = loadKeggPathways()

    # opens up the input and output file.
    f1 = open(input_file, "r")
    f2 = open(output_file, "w")

    # Goes through each line of the predicted protein file.
    for line in f1:

        # Checks to see if the BLAST has a uniprot number associated with it, else moves on to next line
        uniprotID = getUniProtFromBlast(line, threshold)
        if not uniprotID:
            continue
        else:

            # Tries to get the kegg genes associated to uniprot ID, else next line.
            try:
                kegg_genes_list = getKeggGenes(uniprotID)
            except IndexError:
                continue

            # Searches each Kegg entry for a potential KEGG orthology record, else next line.
            for kegg_gene in kegg_genes_list:
                try:
                    kegg_ortho_list = getKeggOrthology(kegg_gene)
                except IndexError:
                    continue
                # Final API call to get associated pathways and maps from KO entry for printing, else continue.
                for kegg_ortho in kegg_ortho_list:
                    kegg_path_list = getKeggPathIDs(kegg_ortho)
                    addKEGGPathways(line, kegg_ortho, kegg_path_list, kegg_paths, f2)

    f1.close()
    f2.close()


def getUniProtFromBlast(blast_line=None, threshold=None):
    """
    Return UniProt ID from the BLAST line if the evalue is below the threshold.
    Returns False if evalue is above threshold.
    """

    # Clean and parse string, and returns false whenever threshold fails or field is incorrect.
    cleaned_line = blast_line.strip()
    blast_fields = cleaned_line.split("\t")
    if float(blast_fields[7]) < float(threshold):
        return blast_fields[1]
    else:
        return False


def loadKeggPathways():
    """
    Return dictionary of key=pathID, value=pathway name from http://rest.kegg.jp/list/pathway/ko
    Example: keggPathways["path:ko00564"] = "Glycerophospholipid metabolism"
    """
    keggPathways = {}
    result = requests.get('https://rest.kegg.jp/list/pathway/ko')
    for entry in result.iter_lines():
        str_entry = entry.decode(result.encoding)  # convert from binary value to plain text
        fields = str_entry.split("\t")
        keggPathways[fields[0]] = fields[1]
    return keggPathways


def getKeggGenes(uniprotID):
    """
    @param uniprotID: uniprot accession number to call from API
    Return a list of KEGG organism:gene pairs for a provided UniProtID.
    """
    kegg_Genes = []
    result = requests.get(f"https://rest.kegg.jp/conv/genes/uniprot:{uniprotID}")
    for entry in result.iter_lines():
        str_entry = entry.decode(result.encoding)
        fields = str_entry.split("\t")
        kegg_Genes.append(fields[1])
    return kegg_Genes


def getKeggOrthology(kegg_Gene):
    """Return Kegg Orthology ID for a KEGG ID."""

    kegg_ortho_genes = []
    result = requests.get(f"https://rest.kegg.jp/link/ko/{kegg_Gene}")
    for entry in result.iter_lines():
        str_entry = entry.decode(result.encoding)  # convert from binary value to plain text
        fields = str_entry.split("\t")
        kegg_ortho_genes.append(fields[1])  # second field is the keggGene value
    return kegg_ortho_genes


def getKeggPathIDs(kegg_Ortho):
    """Return Kegg Path ID for Kegg Orthology ID in the list."""

    kegg_path_genes = []
    result = requests.get(f"https://rest.kegg.jp/link/pathway/{kegg_Ortho}")
    for entry in result.iter_lines():
        str_entry = entry.decode(result.encoding)
        if re.search(r"map", str_entry):
            fields = str_entry.split("\t")
            kegg_path_genes.append(fields[1])
        elif re.search(r"path", str_entry):
            fields = str_entry.split("\t")
            kegg_path_genes.append(fields[1])
        else:
            continue
    return kegg_path_genes


def addKEGGPathways(line=None, kegg_ortho=None, kegg_path_list=None, kegg_paths=None, output_file=None):
    """Print out the output"""
    kegg_path_list_refine =[]
    for entry in kegg_path_list:
        if re.search(r"map", entry):
            continue
        kegg_path_list_refine.append(entry)
    for path in kegg_path_list_refine:
        print(f"{line.strip()}\t{kegg_ortho}\t{path}\t{kegg_paths[path]}", file=output_file)


def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Get the arguments for Assignment 10",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input_file',
                        metavar='INPUT_FILE',  
                        help='need a predicted protein input file',  
                        type=str,  
                        default='/scratch/martinez.jam/module10-jamalemar/data/alignPredicted.txt',
                        required=False
                        )
    # Create a sequential argument (eg. it has to come in the order defined)
    parser.add_argument('-t', '--threshold',  # name of the argument, 
                        metavar='THRESHOLD',  # shorthand to represent the input value
                        help='This is the e-value threshold',  # message to the user, it goes into the help menu
                        type=float,  # type of input expected, could also be int or float
                        default=1e-50,  # default option if no input is given by the user
                        required=False # whether this input must be given by the user, could also be True
                        )
    # Create a flagged argument (eg. input comes after a short "-i" or long "--input" form flag)
    parser.add_argument('-o', '--output_file',
                        # name of the argument, we will later use args.number to get this user input
                        metavar='OUTPUT_FILE',  # shorthand to represent the input value
                        help='File to output results to.',  # message to the user, it goes into the help menu
                        type=str,  # type of input expected, could also be int or float
                        default="/scratch/martinez.jam/module10-jamalemar/data/output.txt",
                        # default option if no input is given by the user
                        required=False  # whether this input must be given by the user, could also be True
                        )

    return parser.parse_args()


if __name__ == "__main__":
    main()

