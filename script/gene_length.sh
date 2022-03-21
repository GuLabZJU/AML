<<eof
awk '{ 
        if(NR==FNR){
            wantedGenes[$1]++ 
        } 
        else{
            if($3=="gene"){
                gsub(/";*/,"",$10); 
                if($10 in wantedGenes){
                    print $10,$5-$4
                } 
            }
        }
      }' genes
eof

awk '{if($3 == "gene"){gsub(/";*/,"",$10); print $10,$5-$4}}' /Users/yinjie/Project/AML/Homo_sapiens.GRCh38.98.gtf > gene.length.tsv

