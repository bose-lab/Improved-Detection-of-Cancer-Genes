random.forest.five<-function(muttab, reps=10000) {
		# 1) Unaffected Residues
		
		pzeroes<-function(gene.table) {
		
			pzero<-function(gene.object){
				#Silent mutations are not used in unaffected residues
				index<-gene.object$Variant_Classification != "Silent" & !is.na(gene.object$Protein_Position)
				if (is.na(gene.object$Protein_Length[1])|sum(index)==0){return(list(PZero=NA))}
				
				#probability of the number of positions with zero mutations, based on poisson distribution for a single zero, and
				#binomial for the number of zeroes. 
				mutnum<-sum(index)
				positions<-trunc(unique(gene.object$Protein_Position[index]))
				plength<-trunc(gene.object$Protein_Length[1])
				
				zeroes<-plength-length(positions)
				
				pzero<-ppois(0, mutnum/plength)
				
				pzeroes<-pbinom(plength-zeroes, plength, 1-pzero)
				return(pzeroes)
				}
			
			return(unlist(as.list(by(gene.table, gene.table$Hugo_Symbol, pzero))))
		}	
		
		# 2) Truncation Rate
		
		loss<-function(gene.table) {
			null<-sum(gene.table$Variant_Classification %in% c("Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"))/length(gene.table$Variant_Classification)
			
			ploss<-function(gene.object, expected=null) {
				observed<-sum(gene.object$Variant_Classification %in% c("Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"))
				p<-pbinom(observed, length(gene.object$Variant_Classification), prob=expected, lower.tail=F)
				return(p)
				}
			
			return(unlist(as.list(by(gene.table, gene.table$Hugo_Symbol, ploss, expected=null))))
			}
		
		# 3) VEST Mean
		
		fis<-function(gene.table, column, iter=10000) {
			levels<-sort(unique(table(gene.table$Hugo_Symbol)))
			null<-list()
			
			#make a multi-tier null for each unique number of mutations
			null[levels]<-lapply(levels, function(x) apply(matrix(sample(gene.table[,column], iter*x, replace=T), ncol=iter), 2, mean, na.rm=T))
			
			#for each gene, how many of the appropriate null are greater than the observed mean?
			return(unlist(as.list(by(gene.table[,column], gene.table$Hugo_Symbol, function(x) sum(null[[length(x)]]>=mean(x, na.rm=T))/iter))))
			}
				
		
		# 4) Patient and Cancer Type Distribution
		
		gof<-function(gene.table, column, iter=10000) {
			levels<-sort(unique(table(gene.table$Hugo_Symbol)))
			null<-list()
			
			gene.table[,column]<-as.numeric(factor(gene.table[,column]))
			expected<-tabulate(gene.table[,column])
			bins<-length(expected)
			
			chisq<-function(obs, ex=expected, nbins=bins) {
				obs<-tabulate(obs, nbins)
				ex<-ex*sum(obs)/sum(ex)
				return(sum((obs-ex)^2/ex))
			}
			
			#make a multi-tier null for each unique level
			null[levels]<-lapply(levels, function(x) apply(matrix(sample(gene.table[,column], iter*x, replace=T), ncol=iter), 2, chisq))
			
			#for each gene, how many of the appropriate null are greater than the observed mean?
			return(unlist(as.list(by(gene.table[,column], gene.table$Hugo_Symbol, function(x) sum(null[[length(x)]]>=chisq(x))/iter))))
		}


		#5) Random Forest 5
		
		rf5<-function(feat.table, label=NULL) {

				#The gold panel, to be used in the absence of user-defined labels.
				ONC<-c("ABL1", "AKT1", "AKT2", "ALK", "BCL6", "BRAF", "CARD11", "CCNE1", "CTNNB1", "EGFR", "ERBB2", "EZH2", "FAS", "FGFR2", 
				"FGFR3", "FLT3", "GNA11", "GNAQ", "GNAS", "HRAS", "IDH1", "JUN", "KDR", "KIT", "KRAS", "MAP2K1", "MAP2K2", "MAP2K4", "MDM2", 
				"MDM4", "MET", "MITF", "MLL", "MYC", "MYCL1", "MYCN", "MYD88", "NFE2L2", "NKX2-1", "NRAS", "PDGFRA", "PIK3CA", "REL", "RET", "RNF43", "SMO", 
				"SOX2", "STAT3", "TERT", "TRAF7", "TSHR")

				TSG<-c("AMER1", "APC", "ATM", "AXIN1", "BAP1", "BRCA1", "BRCA2", "CDH1", "CDKN2A", "CDKN2C", "CEBPA", "CREBBP", "CYLD", 
				"DICER1", "EP300", "FBXW7", "GATA3", "HNF1A", "KDM6A", "MAX", "MEN1", "MLH1", "MSH2", "MSH6", "NF1", "NF2", "NOTCH1", "NOTCH2", 
				"PAX5", "PIK3R1", "PRKAR1A", "PTCH1", "PTEN", "RB1", "SETD2", "SMAD4", "SMARCA4", "SMARCB1", "SOCS1", "STK11", "SUFU", "TET2", 
				"TNFAIP3", "TP53", "TSC1", "TSC2", "VHL", "WT1")
				
				#Returns a feature table with missing or infinite values replaced by the column mean
				impute<-function(gt) {
					imp<-function(x) {
						x[is.infinite(x)]<-mean(x[!is.infinite(x)], na.rm=T)
						x[is.na(x)]<-mean(x, na.rm=T)
						return(x)
						}
					
					return(apply(gt, 2, imp))
					}
				
				feat.table<-impute(feat.table)
				
				#Making a label vector using ONC and TSG if necessary
				if (is.null(label)) {
					label<-rep("UK", dim(feat.table)[1])
					label[rownames(feat.table) %in% ONC]<-"ONC"
					label[rownames(feat.table) %in% TSG]<-"TSG"
					}
				
				
				#60% of minor gene classes are used to train individual trees, with 10x more of the major class.
				frequency<-table(label)
				train.size<-rep(trunc(min(frequency)*0.6), length(frequency))
				train.size[which.max(frequency)]<-train.size[which.min(frequency)]*10

				
				rf<-require("randomForest", quietly=T)	
				if (!rf) {
					print("Without the randomForest package, only individual features can be calculated.")
					return(list(Label=label))
					}
			
				randomforest<-randomForest(feat.table, as.factor(label), ntree=1500, mtry=2, sampsize=train.size)
				
				
				out<-as.data.frame(randomforest$votes)
				colnames(out)<-paste0("Custom_RF_Score_", colnames(out), sep="")
				out$Label<-label
				
				return(out)
				
			}
		
		# 7) Published Forest 5
		
		pf5<-function(feat.table) {
		
			feat.table<-feat.table[,c("Patient.Distribution", "Unaffected.Residues", "VEST.Mean", "Truncation.Rate", "Cancer.Type.Distribution")]

		
			impute<-function(gt) {
				imp<-function(x) {
					x[is.infinite(x)]<-mean(x[!is.infinite(x)], na.rm=T)
					x[is.na(x)]<-mean(x, na.rm=T)
					return(x)
					}
			
			return(apply(gt, 2, imp))
			}
				
			feat.table<-impute(feat.table)
			rf5_model<-NULL
			
			try(load("final_pan-cancer_RF5_model.Rdata"), silent=T)

			if (!is.null(rf5_model)) {
				require(randomForest)
				colnames(feat.table)<-rownames(importance(rf5_model))
				votes<-predict(rf5_model, feat.table, type="prob")
				colnames(votes)<-paste0("Published_RF5_Score_", colnames(votes))
				return(votes)
				}
			else {
				print("final_pan-cancer_RF5_model.Rdata not found. Published model will not be applied.")
				return(NULL)
				}
			
			}


		feat<-list()
		print("Calculating Tests:")
		flush.console()
		
		feat[["Unaffected Residues"]]<-pzeroes(muttab)
		print("1/5 Done")
		flush.console()
		
		feat[["Truncation Rate"]]<-loss(muttab)
		print("2/5 Done")
		flush.console()
		
		feat[["VEST Mean"]]<-fis(muttab, "VEST3", reps)
		print("3/5 Done")
		flush.console()
		
		feat[["Patient Distribution"]]<-gof(muttab, "Sample", reps)
		print("4/5 Done")
		flush.console()
		
		feat[["Cancer Type Distribution"]]<-gof(muttab, "Cancer_Type", reps)
		print("5/5 Done")
		flush.console()

		lab<-list()
		feat<-as.data.frame(feat)
		
		if (any(colnames(muttab)=="Label")) {
			lab<-unique(muttab[,c("Hugo_Symbol", "Label")]) 
			if (any(duplicated(lab$Hugo_Symbol))) {
				print("A gene cannot belong to more than one gene class. Default labels will be used.")
				lab<-list()
				lab$Label<-NULL
				} 
			} 
		


			
		print("Training Custom Model and Applying Published Model")
		flush.console()
		
		custom_pred<-rf5(feat, lab$Label)
		pub_pred<-pf5(feat)
		
		
		if (!is.null(pub_pred)) {
			return(cbind(Hugo_Symbol=unique(muttab$Hugo_Symbol), feat, custom_pred, pub_pred))
			}
		else {
			return(cbind(Hugo_Symbol=unique(muttab$Hugo_Symbol), feat, custom_pred))
			}
		}
	



input<-commandArgs(trailingOnly=T)[1]
		
if (is.na(input)) {
	print("Usage: Rscript model_RF5.r <inputfile>")
	quit()
	} else {
	
	infile<-read.table(input, sep="\t", stringsAsFactors=T, header=T)
	}
	
outfile<-random.forest.five(infile,100)
output<-paste(gsub(".txt", "", input), ".rf5.output.txt", sep="")
write.table(outfile, output, sep="\t", row.names=F, quote=F) 
print(paste("Results file written to ", output, sep=""))
flush.console()
	