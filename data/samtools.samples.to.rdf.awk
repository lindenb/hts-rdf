BEGIN	{
	FS="\t";
	}

	{
	if($1!=".") {
		printf("<foaf:Person rdf:about=\"&u1087;samples/%s\">\n<foaf:name>%s</foaf:name>\n",$2,$2);
		printf("</foaf:Person>\n");
		}
	printf("<u:Bam rdf:about=\"file://%s\"><u:filename>%s</u:filename>\n",$3,$3);
	if($1!=".") printf("<u:sample rdf:resource=\"&u1087;samples/%s\"/>\n",$2);
	if($3!="." && $3!="") printf("<u:reference rdf:resource=\"&u1087;references/%s\"/>",$4);
	printf("</u:Bam>\n");
	}
