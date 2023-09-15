#
# This makefile is used to generate README.md ,  knowledge.rdf and knowledge.svg
#

SHELL=/bin/bash
#JENA=${HOME}/package/apache-jena-4.8.0
JENA=${HOME}/packages/apache-jena-4.8.0/
JENA_BIN_DIR=$(JENA)/bin
ARQ=$(JENA_BIN_DIR)/arq
JSANDBOX=../jsandbox

all: README.md


README.md : data/README.m4 \
		knowledge.rdf \
		knowledge.svg \
		data/diseases.svg \
		data/species.svg \
		TMP/species.01.out \
		TMP/diseases.01.out \
		TMP/samples.01.out \
		TMP/bams.01.out \
		TMP/vcfs.01.out \
		TMP/references.rdf \
		data/tsv2markdown.awk
	m4 -P < $< > $(addsuffix .tmp,$@)
	mv $(addsuffix .tmp,$@) $@



knowledge.svg : knowledge.rdf
	mkdir -p $(dir $@)
	java -jar $(JSANDBOX)/dist/rdf2graph.jar $<  | dot -T svg -o $@

data/diseases.svg : data/diseases.rdf
	mkdir -p $(dir $@)
	java -jar $(JSANDBOX)/dist/rdf2graph.jar $<  | dot -T svg -o $@

data/species.svg : data/species.rdf
	mkdir -p $(dir $@)
	java -jar $(JSANDBOX)/dist/rdf2graph.jar $<  | dot -T svg -o $@



TMP/species.01.out: knowledge.rdf data/query.species.01.sparql
	mkdir -p $(dir $@)
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)  --results=TSV > $@

TMP/diseases.01.out: knowledge.rdf data/query.diseases.01.sparql
	mkdir -p $(dir $@)
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)  --results=TSV > $@

TMP/samples.01.out: knowledge.rdf data/query.samples.01.sparql
	mkdir -p $(dir $@)
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)  --results=TSV > $@

TMP/bams.01.out: knowledge.rdf data/query.bams.01.sparql
	mkdir -p $(dir $@)
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)  --results=TSV > $@

TMP/vcfs.01.out: knowledge.rdf data/query.vcfs.01.sparql
	mkdir -p $(dir $@)
	$(ARQ) --data=$(word 1,$^) --query=$(word 2,$^)  --results=TSV > $@


knowledge.rdf : TMP/references.rdf data/species.rdf TMP/bams.rdf TMP/bams.rdf data/diseases.rdf data/samples.rdf TMP/vcf2ref.rdf  TMP/vcf2samples.rdf
	$(JENA_BIN_DIR)/riot --formatted=RDFXML $^ > $@



#
# convert REF to rdf
#
TMP/references.rdf : data/references.tsv data/header.rdf.part data/footer.rdf.part
	mkdir -p $(dir $@)
	cat data/header.rdf.part > $@
	tail -n+2 $< |\
		awk -F '\t' '{printf("<u:Reference rdf:about=\"&u1087;references/%s\">\n\t<u:genomeId>%s</u:genomeId>\n\t<u:filename>%s</u:filename>\n",$$2,$$2,$$1);if($$4!="") printf("\t<u:taxon rdf:resource=\"http://purl.uniprot.org/taxonomy/%s\"/>\n",$$4); printf("</u:Reference>\n");}'  >> $@
	cat data/footer.rdf.part >> $@
#
# convert samples+REF to RDF
#
TMP/bams.rdf :  TMP/bams.txt data/header.rdf.part data/footer.rdf.part data/samtools.samples.to.rdf.awk
	mkdir -p $(dir $@)
	cat data/header.rdf.part > $@
	awk -F '\t' -f data/samtools.samples.to.rdf.awk $< >> $@
	cat data/footer.rdf.part >> $@



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


TMP/vcf2ref.rdf :TMP/vcfs.md5.txt TMP/references.md5 data/header.rdf.part data/footer.rdf.part
	mkdir -p $(dir $@)
	cat data/header.rdf.part > $@
	join -t $$'\t' -1 2 -2 5 $(word 1,$^) $(word 2,$^) |\
		awk -F '\t' '{printf("<u:Vcf rdf:about=\"file://%s\"><u:filename>%s</u:filename><u:reference rdf:resource=\"&u1087;references/%s\"/></u:Vcf>",$$2,$$2,$$4); }' >> $@
	cat data/footer.rdf.part >> $@

TMP/vcf2samples.rdf: TMP/vcfs.txt  data/header.rdf.part data/footer.rdf.part
	mkdir -p $(dir $@)
	cat data/header.rdf.part > $@
	cat $< | while read F; do bcftools query -l "$${F}" | awk -vVCF="$$F" 'BEGIN {printf("<u:Vcf rdf:about=\"file://%s\"><u:filename>%s</u:filename>",VCF,VCF); } {printf("<u:sample rdf:resource=\"&u1087;samples/%s\"/>",$$1);} END {printf("</u:Vcf>");}' >> $@ ; done
	cat data/footer.rdf.part >> $@

# extract contig/length from vcf header, calc the md5
TMP/vcfs.md5.txt: TMP/vcfs.txt
	mkdir -p $(dir $@)
	(cat $<| while read V ; \
			do echo -en "$${V}\t" && \
			bcftools view --header-only "$${V}"  | awk  -F '[=,<>]' '/^##contig/ {printf("%s\t%s\n",$$4,$$6);}'   | md5sum | cut -d ' ' -f1 ; done) | sort -t $$'\t' -k2,2 > $@


## find all the VCF in data
TMP/vcfs.txt:
	mkdir -p $(dir $@)
	find data -type f \( -name "*.vcf.gz" -o -name "*.bcf" -o -name "*.vcf" \) | sort > $@

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

TMP/references.md5: TMP/sorted.refs.txt
	mkdir -p $(dir $@)
	cut -f 1 $< | while read FA; do echo -ne "$${FA}\t" && cut -f 1,2 "$${FA}.fai" | md5sum | cut -d ' ' -f 1; done > $(addsuffix .tmp.a,$@)
	join -t $$'\t' -1 1 -2 1  $< $(addsuffix .tmp.a,$@) | sort -t $$'\t' -k5,5 > $@
	rm -f  $(addsuffix .tmp.a,$@)

clean:
	rm -rf TMP
