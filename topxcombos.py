#!/usr/bin/python

#=============================================================================
# Top X Combos
# Retrieves the top X percent of scored overhang combinations from an output
#   file from GGSOFT.
#
# Copyright 2013 Michael Ting
# https://github.com/michaelting
#=============================================================================

from argparse import ArgumentParser

def total_lines(infile):
    """
    Calculates the total number of lines in an input file
    
    Input:
        infile      - A file containing the output from ggsoft.py
    Output:
        num_lines   - The total number of lines in infile
    """
    template = open(infile)
    num_lines = sum(1 for line in template)
    template.close()
    
    return num_lines
    
def lines_to_return(percent, total_lines):
    """
    Calculates the number of lines to write to the output file based on the
    percentage of combinations specified by the user
    
    Input:
        percent     - The top X percent of combinations to be written to 
                      an output file.
        total_lines - Total number of lines (combinations) in the input file
    Output:
        lines_returned  - Number of combinations returned based on percentage
    """
    decimal_val = 0.01*percent
    lines_returned = int(decimal_val * total_lines)
    return lines_returned
    
def write_output(infile, outfile, percent, lines_to_return, total_lines):
    """
    Creates the output file with the top X percent of combinations from a
    ggsoft.py output file
    
    Input:
        infile          - The input file corresponding to results from ggsoft.py
        outfile         - The file where the top X percent of combinations is
                          written to
        percent         - The top X percent specified by the user
        lines_to_return - Number of lines returned specified by the percent
        total_lines     - Total number of lines in the input file
    Output:
        newfile         - A file containing the top X percent of combinations
                          from the input ggsoft.py
    """
    
    newfile = open(outfile, "w")
    
    newfile.write("Top %f percent of combinations from %s\n" % (percent, infile))
    newfile.write("%d of %d combinations selected\n" % (lines_to_return, total_lines))

    template = open(infile)
    
    lineindex = 0    
    for line in template:
        if lineindex == lines_to_return:
            break
        newfile.write(line)
        lineindex += 1
        
    template.close()
    newfile.close()

def main():
    
    # read command-line arguments
    parser = ArgumentParser(description="Indicate top X percent of overhang combinations to retrieve")
    
    parser.add_argument("infile", help="Input file with results from GGSOFT")
    parser.add_argument("outfile", help="Output file with top X percent of overhang combinations")
    parser.add_argument("percent", help="Top X percent of combinations to retrieve")
    
    args = parser.parse_args()
    
    infile = args.infile
    outfile = args.outfile
    percent = float(args.percent)
    
    # find total number of lines in the input file
    num_lines = total_lines(infile)
    
    # calculate number of lines to copy to new file based on percent specified
    lines_returned = lines_to_return(percent, num_lines)

    # create output file with top X percent of combinations
    write_output(infile, outfile, percent, lines_returned, num_lines)

if __name__ == "__main__":
    main()