BEGIN	{
	FS="\t";
	}

	{
	for(i=1;i<=NF;i++) {
		S = $i;
		if(NR==1) gsub(/\?/,"",S);
		printf("| %s ",S);
		}
	printf(" |\n");
	if(NR==1) {
		for(i=1;i<=NF;i++) {
			printf("|-----");
			}
		printf("|\n");
		}
	}
