#!/usr/bin/env python

import sys
import argparse
import os
import textwrap

def version():
	v = "----"
	return v

def warning(*objs):
	print("WARNING: ", *objs, end='\n', file=sys.stderr)
	sys.exit()

def get_parser():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=textwrap.dedent(version())
	)

	parser.add_argument('-i','--input', help='input gene info', type=str)
	parser.add_argument('-o', '--output', help='output folder', type=str)
	parser.add_argument('-p', '--prefix', help='output prefix', type=str)

	parser.add_argument('-f', '--fasta', help='fasta folder', type=str)
	parser.add_argument('-g', '--gff', help='gff folder', type=str)
	parser.add_argument('-u', '--up', help='gene upstream (default=10000)', type=int, default=10000)
	parser.add_argument('-d', '--down', help='gene downstream (default=5000)', type=int, default=5000)

	parser.add_argument('-a', '--all', help='run one step (default=True)', choices=[True, False], default=True)

	return parser

def mapping1(info, fasta, gff, up, down, o ,name, onestep):
	os.system("mkdir -p " + o)
	with open (info) as f:
	    for i in f:
	        i = i[:-1].split()
	        gene = i[0]
	        g = i[2]

	        fa = fasta + "/" + g + "*.fa"
	        gff3 = gff + "/" + g + "*" + i[1] + "*gff3.gz"
	        out = os.popen("ls " + fa).read()[:-1].split("\n")
	        print (out)
	        if len(out) == 1 and not out[0].startswith("ls"):
	            fa = out[0]
	        else:
	            print (str(i) + " error!!!!!!!!!")
	        out = os.popen("ls " + gff3).read()[:-1].split("\n")
	        if len(out) == 1 and not out[0].startswith("ls"):
	            gff3 = out[0]
	        else:
	            print (str(i) + " error!!!!!!!!!")

	        
	        gene_info = os.popen("""zcat """ + gff3 + ' |grep "' + gene + '"' + """|awk '{if($3=="gene"){print $1"\\t"$4"\\t"$5"\\t"$7}}'""").read()[:-1].split("\t")
	        cmd = ("""zcat """ + gff3 + ' |grep "' + gene + '"' + """ > """ + o + "/" + gene + ".gff")
	        print (cmd)
	        os.system(cmd)

	        chro = gene_info[0]
	        if gene_info[3] == "+":
	            start = gene_info[1]
	            end = gene_info[2]
	            s = int(start) - int(up)
	            e = int(end) + int(down)
	            os.system("awk 'BEGIN{print " + '"' + chro + '\\t' + str(s) + '\\t' + str(e) + '\\t' + gene + '_' + g + '"' + ";exit}{}'|bedtools getfasta -fi " + fa + " -bed - -name+ -fo " + o + "/" + gene + ".fa")
	        else:
	            start = gene_info[2]
	            end = gene_info[1]
	            s = int(start) + int(up)
	            e = int(end) - int(down)
	            os.system("awk 'BEGIN{print " + '"' + chro + '\\t' + str(e) + '\\t' + str(s) + '\\t' + gene + '_' + g + '"' + ";exit}{}'|bedtools getfasta -fi " + fa + " -bed - -name+ -fo " + o + "/" + gene + ".fa")

	os.system("cat " + o + "/*fa > " + name + ".fa")
	os.system("muscle -in " + name + ".fa -out " + name + ".2.fa -tree2 " + name + ".2.nwk -maxiters 2")

	if onestep:
		os.system("Rscript plot_syn.r --nwk " + name + ".2.nwk --gffdir " + name + "/ --up " + str(up) + " --down " + str(down) + " " + name + ".2.fa " + name)

def main():
	parser = get_parser()
	args = vars(parser.parse_args())

	if args['input'] is not None:
		print(version())

	
	mapping1(args['input'], args['fasta'], args['gff'], args['up'], args['down'], args['output'], args['prefix'], args['all'])


if __name__ == "__main__":
	main()
