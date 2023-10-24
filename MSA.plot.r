args<-commandArgs(T)
library(pheatmap)
data=read.table(args[1],header=0,row.names=1)
data_matrix=as.matrix(data)
data_mark=as.matrix(data)
for(i in 1:nrow(data_matrix)){
	for(j in 1:ncol(data_matrix)){
		if(data_matrix[i,j] == 1){
			data_mark[i,j] <- "A";
		}
		else if(data_matrix[i,j] == 2){
			data_mark[i,j] <- "T";
		}
		else if(data_matrix[i,j] == 3){
			data_mark[i,j] <- "C";
		}
		else if(data_matrix[i,j] == 4){
			data_mark[i,j] <- "G";
		}
		else if(data_matrix[i,j] == 0){
			data_mark[i,j] <- "-" ;
		}
	}
}

pheatmap(data,
	 cluster_col = FALSE,
	 cluster_row = FALSE,
	 legend = FALSE,
	 show_rownames=T,
	 show_colnames=F,
	 cellheight=10,
	 cellwidth=10,
	 filename =args[2],
	 display_numbers=data_mark
)
