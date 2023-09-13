m4_define(`md_pre', m4_changequote([,])[m4_changequote([,])```[$1]m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`md_code', m4_changequote([,])[m4_changequote([,])`[$1]`m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`anchor',`[$1]($2)')m4_dnl
m4_define(`localfile',`anchor($1,$1)')m4_dnl
# hts-rdf

Author: Pierre Lindenbaum PhD.


Here are a few notes about Managing  sequencing data with RDF. I want to keep track of the samples, BAMs, references, diseases etc..  used in my lab.

  - I don't want to use a sql database.
  - I don't want to join too many tab delimited files.
  - I want to use a controlled vocabulary to define things like a disease, etc...
  - In this document I won't explain what are md_code(RDF) or md_code(SPARQL).
  - I use the md_code(RDF+XML) notation because I 'm used to work with md_code(XML).
  - I create a namespace for my lab:  md_code(https://umr1087.univ-nantes.fr/rdf/) and a md_code(XML) entity for this namespace: md_code(&u1087;).
  - I tried to reuse existing ontologies (e.g. md_code(foaf:Person) for samples) as much as I can, but sometimes I created my own classes and properties.


# Building the RDF GRAPH

## Species

I manually wrote localfile(data/species.rdf) defining the species used in my lab.
We will use md_code(rdf:subClassOf) to find organisms that are a sub-species of a taxon in the NCBI taxonomy tree.

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  data/species.rdf)m4_dnl
md_pre


## Diseases

I manually wrote localfile(data/diseases.rdf) defining the diseases used in my lab.
We will use md_code(rdf:subClassOf) to find organisms that are a sub-disease in a disease ontology tree.

md_pre(rdf)
(...)
m4_syscmd(tail -n+24  data/diseases.rdf)m4_dnl
md_pre


## References / FASTA / Genomes

I manually wrote localfile(data/references.tsv) a tab delimited text file defining each FASTA reference genome available on my cluster. The taxon id will be used to retrive the species associated to a FASTA file.

md_pre
m4_syscmd(cat data/references.tsv | column -t -s '	')m4_dnl
md_pre

using awk, the table is transformed into **RDF**:

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

## BAM files

BAM file contains the sample names in their read-groups; We use md_code(samtools samples) to extract the samples, the reference and the path of each BAM file.
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

## Combing all the RDF chunks

[jena/rio](https://jena.apache.org/) is used to merge md_code(RDF) files into localfile(knowledge.rdf)

md_pre(bash)
riot --formatted=RDFXML TMP/references.rdf data/species.rdf TMP/bams.rdf data/diseases.rdf data/samples.rdf > knowledge.rdf
md_pre


# Querying the GRAPH

[jena/arq](https://jena.apache.org/) is used to run the md_code(SPARQL) queries.

md_pre(bash)
arq --data=knowledge.rdf --query=querysparql
md_pre

## Example

> show me the species that are a sub-taxon of "Homo"

query localfile(data/query.species.01.sparql):

md_pre(sparql)
(...)
m4_syscmd(grep -v ^PREFIX  data/query.species.01.sparql)m4_dnl
md_pre

output:

md_pre
m4_syscmd(cat  TMP/species.01.out)m4_dnl
md_pre

## Example

> show the diseases that are a sub disease of **COVID-19**.

query localfile(data/query.diseases.01.sparql):

md_pre(sparql)
(...)
m4_syscmd(grep -v ^PREFIX  data/query.diseases.01.sparql)m4_dnl
md_pre

output:

md_pre
m4_syscmd(cat  TMP/diseases.01.out)m4_dnl
md_pre


## Example

> find the samples , their children, parents , diseases

query localfile(data/query.samples.01.sparql) :

md_pre(sparql)
(...)
m4_syscmd(grep -v ^PREFIX  data/query.samples.01.sparql)m4_dnl
md_pre

output:

md_pre
m4_syscmd(cat  TMP/samples.01.out)m4_dnl
md_pre


## Example

> find the bam , their reference, samples , etc.. 

query localfile(data/query.bams.01.sparql) :

md_pre(sparql)
(...)
m4_syscmd(grep -v ^PREFIX  data/query.bams.01.sparql)m4_dnl
md_pre

output:

md_pre
m4_syscmd(cat  TMP/bams.01.out)m4_dnl
md_pre


# The Graph

and here is the RDF graph

![knowledge.svg](knowledge.svg)
