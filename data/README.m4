m4_define(`md_pre', m4_changequote([,])[m4_changequote([,])```[$1]m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`md_code', m4_changequote([,])[m4_changequote([,])`[$1]`m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`anchor',`[$1]($2)')m4_dnl
m4_define(`localfile',`anchor($1,$1)')m4_dnl
m4_define(`sparql_example',`

query localfile($1) :

md_pre(sparql)
(...)
m4_syscmd(grep -v ^PREFIX  $1)m4_dnl
md_pre

execute:

md_pre(bash)
arq --data=knowledge.rdf --query=$1 > $2
md_pre

output localfile($2):

md_pre
m4_syscmd(cat $2)m4_dnl
md_pre


')m4_dnl
# hts-rdf

Author: Pierre Lindenbaum PhD.


Here are a few notes about Managing  sequencing data with RDF. I want to keep track of the samples, BAMs, references, diseases etc..  used in my lab.

  - This document is auto-generated using a localfile(Makefile). Do not edit it.
  - I don't want to use a md_code(SQL) database.
  - I don't want to join too many tab delimited files.
  - I want to use a controlled vocabulary to define things like diseases, organims, etc...
  - This document is NOT a tutorial for md_code(RDF) or md_code(SPARQL).
  - I use the md_code(RDF+XML) notation because I 'm used to work with md_code(XML).
  - I created a namespace for my lab:  md_code(https://umr1087.univ-nantes.fr/rdf/) and a md_code(XML) entity for this namespace: md_code(&u1087;).
  - I tried to reuse existing ontologies (e.g. md_code(foaf:Person) for samples) as much as I can, but sometimes I created my own classes and properties.
  - I'm not an expert of md_code(SPARQL) or md_code(RDF)
  - Required tools are (jena)[https://jena.apache.org/download/], md_code(bcftools= (for md_code(VCFs)), md_code(samtools) (for md_code(BAMs)), md_code(awk).


# Building the RDF GRAPH

## Species

I manually wrote localfile(data/species.rdf) defining the species used in my lab.
We will use md_code(rdf:subClassOf) to find organisms that are a sub-species of a taxon in the NCBI taxonomy tree.

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  data/species.rdf)m4_dnl
md_pre


## Diseases / Phenotypes

I manually wrote localfile(data/diseases.rdf) defining the diseases used in my lab.
We will use md_code(rdf:subClassOf) to find diseases that are a sub-disease in a disease ontology tree.

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  data/diseases.rdf)m4_dnl
md_pre


## References / FASTA / Genomes

I manually wrote localfile(data/references.tsv) a tab delimited text file defining each md_code(FASTA) reference genome available on my cluster.
The taxon id will be used to retrive the species associated to a md_code(FASTA) file.

md_pre
m4_syscmd(cat data/references.tsv | column -t -s '	')m4_dnl
md_pre

The table is transformed into md_code(RDF) using md_code(awk):

md_pre(bash)
tail -n+2 data/references.tsv |\
	awk -F '\t' '{printf("<u:Reference rdf:about=\"&u1087;references/%s\">\n\t<u:genomeId>%s</u:genomeId>\n\t<u:filename>%s</u:filename>\n",$2,$2,$1);if($4!="") printf("\t<u:taxon rdf:resource=\"http://purl.uniprot.org/taxonomy/%s\"/>\n",$4); printf("</u:Reference>\n");}'
md_pre

output:

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  TMP/references.rdf)m4_dnl
md_pre

## Samples


I manually wrote localfile(data/samples.rdf) defining the samples sequenced in my lab.
This is where we can define the gender, associate a sample to a diseases and where we can define the familial relations.
The Class md_code(foaf:Group) is used to create a group of samples.

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  data/samples.rdf)m4_dnl
md_pre

## VCF BCF

### VCF and genomes

for md_code(VCF) files, we need to associate a md_code(VCF) and the reference genome:
md_code(Chromosome) and md_code(length) are extracted from the references, we calculate the  md_code(md5) checksum and we sort on  md_code(md5).

md_pre(bash)
tail -n+2 data/references.tsv | sort -T TMP -t $'\t' -k1,1 > TMP/sorted.refs.txt
cut -f 1 TMP/sorted.refs.txt | while read FA; do echo -ne "${FA}\t" && cut -f 1,2 "${FA}.fai" | md5sum | cut -d ' ' -f 1; done > TMP/references.md5.tmp.a
join -t $'\t' -1 1 -2 1  TMP/sorted.refs.txt TMP/references.md5.tmp.a | sort -t $'\t' -k5,5 > TMP/references.md5
rm -f  TMP/references.md5.tmp.a
md_pre

For each md_code(VCF), the header is extracted, we extract the the md_code(chromosome) and md_code(length) of the md_code(contig) lines, we calculate the md_code(md5) checksum and we sort on  md_code(md5).

md_pre(bash)
find data -type f \( -name "*.vcf.gz" -o -name "*.bcf" -o -name "*.vcf" \) | sort > TMP/vcfs.txt
(cat TMP/vcfs.txt| while read V ; \
		do echo -en "${V}\t" && \
		bcftools view --header-only "${V}"  | awk  -F '[=,<>]' '/^##contig/ {printf("%s\t%s\n",$4,$6);}'   | md5sum | cut -d ' ' -f1 ; done) | sort -t $'\t' -k2,2 > TMP/vcfs.md5.txt
md_pre

we join both files on md_code(md5) and we convert to md_code(RDF) using md_code(awk):

md_pre(bash)
cat data/header.rdf.part > TMP/vcf2ref.rdf
join -t $'\t' -1 2 -2 5 TMP/vcfs.md5.txt TMP/references.md5 |\
	awk -F '\t' '{printf("<u:Vcf rdf:about=\"file://%s\"><u:filename>%s</u:filename><u:reference rdf:resource=\"&u1087;references/%s\"/></u:Vcf>",$2,$2,$4); }' >> TMP/vcf2ref.rdf
cat data/footer.rdf.part >> TMP/vcf2ref.rdf
md_pre

### VCF and samples

to link the  md_code(VCF) files and the  sample, we use  md_code(bcftools query -l) to extract the samples and we convert to md_code(RDF) using md_code(awk):

md_pre(bash)
find data -type f \( -name "*.vcf.gz" -o -name "*.bcf" -o -name "*.vcf" \) | sort > TMP/vcfs.txt
cat data/header.rdf.part > TMP/vcf2samples.rdf
cat TMP/vcfs.txt | while read F; do bcftools query -l "${F}" | awk -vVCF="$F" 'BEGIN {printf("<u:Vcf rdf:about=\"file://%s\"><u:filename>%s</u:filename>",VCF,VCF); } {printf("<u:sample rdf:resource=\"&u1087;samples/%s\"/>",$1);} END {printf("</u:Vcf>");}' >> TMP/vcf2samples.rdf ; done
cat data/footer.rdf.part >> TMP/vcf2samples.rdf
md_pre


## BAM files

md_code(BAM) file contains the sample names in their read-groups; We use md_code(samtools samples) to extract the samples, the reference and the path of each md_code(BAM) file.
localfile(data/samtools.samples.to.rdf.awk) is used to convert the output of md_code(samtools samples)  to md_code(RDF).


md_pre(bash)
find ${PWD}/data -type f -name "*.bam" |\
	samtools samples -F TMP/references.txt |\
	sort -T TMP -t $'\t' -k3,3 |\
	join -t $'\t' -1 3 -2 1 - TMP/sorted.refs.txt > TMP/bams.txt

cat data/header.rdf.part > TMP/bams.rdf

awk -F '\t' -f data/samtools.samples.to.rdf.awk TMP/bams.txt >> TMP/bams.rdf

cat data/footer.rdf.part >> TMP/bams.rdf
md_pre

the output:

md_pre(rdf)
(...)
m4_syscmd(xmllint --format TMP/bams.rdf | grep S5 -A 2)m4_dnl
(...)
md_pre

## Combining all the RDF chunks

[jena/rio](https://jena.apache.org/) is used to merge md_code(RDF) files into localfile(knowledge.rdf)

md_pre(bash)
riot --formatted=RDFXML TMP/references.rdf data/species.rdf TMP/bams.rdf data/diseases.rdf data/samples.rdf TMP/vcf2ref.rdf TMP/vcf2samples.rdf > knowledge.rdf
md_pre


# Querying the GRAPH

[jena/arq](https://jena.apache.org/) is used to run the md_code(SPARQL) queries.

md_pre(bash)
arq --data=knowledge.rdf --query=querysparql
md_pre

## Example

> show me the species that are a sub-taxon of "Homo"

sparql_example(data/query.species.01.sparql,TMP/species.01.out)

## Example

> show the diseases that are a sub disease of **COVID-19**.

sparql_example(data/query.diseases.01.sparql,TMP/diseases.01.out)

## Example

> find the samples , their children, parents , diseases

sparql_example(data/query.samples.01.sparql,TMP/samples.01.out)

## Example

> List all the md_code(VCF) files and their samples, at least containing the sample "S1"

sparql_example(data/query.vcfs.01.sparql,TMP/vcfs.01.out)

## Example

> find the bam , their reference, samples , etc..

sparql_example(data/query.bams.01.sparql,TMP/bams.01.out)

# The Graph

and here is the md_code(RDF) graph as a md_code(SVG) document:

![knowledge.svg](knowledge.svg)
