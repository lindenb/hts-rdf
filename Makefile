SHELL=/bin/bash
JENA=${HOME}/package/apache-jena-4.8.0
JENA_BIN_DIR=$(JENA)/bin
ARQ=$(JENA_BIN_DIR)/arq


test:data/species.rdf data/query.species.01.sparql
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)

#
# convert REF to rdf
#
TMP/references.rdf : data/references.tsv data/header.rdf.part data/footer.rdf.part
	mkdir -p $(dir $@)
	cat data/header.rdf.part > $@
	cat data/header.rdf.part > $@
	tail -n+2 $< |\
		awk -F '\t' '{printf("<u:Reference rdf:about=\"&u1087;references/%s\">\n\t<u:genomeId>%s</u:genomeId>\n\t<u:filename>%s</u:filename>\n",$$2,$$2,$$1);if($$4!="") printf("\t<u:taxon rdf:resource=\"http://purl.uniprot.org/taxonomy/%s\"/>\n",$$4); printf("</u:Reference>\n");}'  >> $@
	cat data/footer.rdf.part >> $@
#
# convert samples+REF to RDF
#
TMP/bams.rdf :  TMP/bams.txt data/header.rdf.part data/footer.rdf.part

#
# build table linking bam and samples
#
TMP/bams.txt : TMP/references.txt TMP/sorted.refs.txt
	mkdir -p $(dir $@)
	# find all bam
	# extract sample name and REF using samtools samples
	# sort on REF
	find ${PWD}/data -type f -name "*.bam" |\
		samtools samples -F $< |\
		sort -T TMP -t $$'\t' -k3,3 |\
		join -t $$'\t' -1 3 -2 1 - TMP/sorted.refs.txt > $@

#
# references for samtools samples
#
TMP/references.txt : data/references.tsv
	mkdir -p $(dir $@)
	tail -n+2 $< | cut -f 1 > $@

#
# references ordered by FASTA
#
TMP/sorted.refs.txt:  data/references.tsv
	mkdir -p $(dir $@)
	tail -n+2 $< | sort -T TMP -t $$'\t' -k1,1 > $@
