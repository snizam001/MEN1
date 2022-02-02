Download data: PRJNA772915. Go through the article to download published datasets. 

#-- quality control filtering
#__________________________________
mkdir Trim_galore

for file in MLL1 MEN1 JunD PolII-RPB1 MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K
do

	trim_galore --paired $file-Exp_1.fastq.gz $file-Exp_2.fastq.gz --output_dir Trim_galore/$file --j 1  &
	trim_galore --paired $file-Ctrl_1.fastq.gz $file-Ctrl_2.fastq.gz --output_dir Trim_galore/$file --j 1 &

done 

#-- alignment 
#__________________________________
mkdir Bamfiles
mkdir SpikeIn
human="/home/sheikh/Databases/hg38/GenCode/GRCh38.p13"
spike="/home/sheikh/Databases/Drosophila_genome/Drosophila_melanogaster"

for file in MLL1 MEN1 JunD PolII-RPB1 MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K
do
	bowtie2 -x $human --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 -p 30 -1 Trim_galore/$file/$file-Exp_1_val_1.fq.gz -2 Trim_galore/$file/$file-Exp_2_val_2.fq.gz | samtools sort -@ 20 -O BAM -o Bamfiles/$file\_expr.bam - 
	samtools index Bamfiles/$file\_expr.bam

	bowtie2 -x $human --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 -p 30 -1 Trim_galore/$file/$file-Ctrl_1_val_1.fq.gz -2 Trim_galore/$file/$file-Ctrl_2_val_2.fq.gz | samtools sort -@ 20 -O BAM -o Bamfiles/$file\_ctrl.bam - 
	samtools index Bamfiles/$file\_ctrl.bam

	bowtie2 -x $spike --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 -p 30 -1 Trim_galore/$file/$file-Exp_1_val_1.fq.gz -2 Trim_galore/$file/$file-Exp_2_val_2.fq.gz | samtools sort -@ 20 -O BAM -o SpikeIn/$file\_expr.bam - 
	samtools index SpikeIn/$file\_expr.bam

	bowtie2 -x $spike --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 -p 30 -1 Trim_galore/$file/$file-Ctrl_1_val_1.fq.gz -2 Trim_galore/$file/$file-Ctrl_2_val_2.fq.gz | samtools sort -@ 20 -O BAM -o SpikeIn/$file\_ctrl.bam - 
	samtools index SpikeIn/$file\_ctrl.bam
	
	samtools flagstat Bamfiles/$file\_expr.bam > Bamfiles/$file\_expr.flagstat.txt 
	samtools flagstat Bamfiles/$file\_ctrl.bam > Bamfiles/$file\_ctrl.flagstat.txt 
	samtools flagstat SpikeIn/$file\_expr.bam > SpikeIn/$file\_expr.flagstat.txt 
	samtools flagstat SpikeIn/$file\_ctrl.bam > SpikeIn/$file\_ctrl.flagstat.txt 
	
done 

#-- calculation of the coverage
#__________________________________________
mkdir bamcoverage

for file in MLL1 JunD PolII-RPB1 
do

bamCoverage -b Bamfiles/$file\_expr.bam -o ./bamcoverage/$file\_expr.bw -bl /home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed -p 30 --effectiveGenomeSize 2913022398

done 

#-- Selection of equal number of reads and coverage
#__________________________________________________
#- 
min=$(grep read1 Bamfiles/MEN1*txt  | grep _expr | awk -F ":" '{print $2}'  | awk '{print $1}'  | sort -nk1  | head -n1 )

for file in MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K MEN1
do
	c=$(grep read1 Bamfiles/$file\_expr.flagstat.txt | awk '{print $1}' )
	ratio=$(echo "$min/$c" | bc -l)
	echo $ratio
	sambamba view -h -t 30 -s $ratio -f bam --subsampling-seed=786 Bamfiles/$file\_expr.bam -o Bamfiles/$file\_expr.selected.bam
	samtools flagstat Bamfiles/$file\_expr.selected.bam  > Bamfiles/$file\_expr.selected.flagstat.txt
	sambamba view -h -t 30 -s $ratio -f bam --subsampling-seed=786 Bamfiles/$file\_ctrl.bam -o Bamfiles/$file\_ctrl.selected.bam 
	samtools flagstat Bamfiles/$file\_ctrl.selected.bam  > Bamfiles/$file\_ctrl.selected.flagstat.txt
	bamCoverage -b Bamfiles/$file\_expr.selected.bam -o ./bamcoverage/$file\_expr.bw -bl /home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed -p 30 --effectiveGenomeSize 2913022398
	
done 


#--Correlation plot
#__________________________________
mkdir CorrelationPlot

multiBamSummary bins -b Bamfiles/MEN1_expr.selected.bam Bamfiles/MEN1-E408Q_expr.selected.bam Bamfiles/MEN1-E255K_expr.selected.bam Bamfiles/MEN1-E359K_expr.selected.bam Bamfiles/MEN1-R52G_expr.selected.bam -o ./CorrelationPlot/out.gz -l "MEN1(WT)" E408Q E255K E359K R52G -bl ~/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed -p 35
plotCorrelation -in ./CorrelationPlot/out.gz -c pearson -p heatmap -o ./CorrelationPlot/out.jpeg --colorMap RdYlBu --plotNumbers  --removeOutliers --plotHeight 6 --plotWidth 7

#-- peak Calling
#__________________________________
mkdir TagDir

#- (preparation: for the experiment)
for file in MLL1 JunD PolII-RPB1
do

	makeTagDirectory TagDir/$file\_expr Bamfiles/$file\_expr.bam

done

for file in MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K
do

	makeTagDirectory TagDir/$file\_expr Bamfiles/$file\_expr.selected.bam

done

#- (preparation: for the control)
for file in MLL1 JunD PolII-RPB1
do
	h_ex=$(fgrep read1 Bamfiles/$file\_expr.flagstat.txt | cut -d " " -f1)
	h_ct=$(fgrep read1 Bamfiles/$file\_ctrl.flagstat.txt | cut -d " " -f1)
	s_ex=$(fgrep read1 SpikeIn/$file\_expr.flagstat.txt  | cut -d " " -f1)
	s_ct=$(fgrep read1 SpikeIn/$file\_ctrl.flagstat.txt  | cut -d " " -f1)
	norm=$(echo "$h_ct*(($s_ct/$h_ct)/($s_ex/$h_ex))" | bc -l)
	makeTagDirectory TagDir/$file\_ctrl Bamfiles/$file\_ctrl.bam -totalReads $norm

done

for file in MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K
do
	h_ex=$(fgrep read1 Bamfiles/$file\_expr.selected.flagstat.txt | cut -d " " -f1)
	h_ct=$(fgrep read1 Bamfiles/$file\_ctrl.selected.flagstat.txt | cut -d " " -f1)
	s_ex=$(fgrep read1 SpikeIn/$file\_expr.flagstat.txt  | cut -d " " -f1)
	s_ct=$(fgrep read1 SpikeIn/$file\_ctrl.flagstat.txt  | cut -d " " -f1)
	norm=$(echo "$h_ct*(($s_ct/$h_ct)/($s_ex/$h_ex))" | bc -l)
	makeTagDirectory TagDir/$file\_ctrl Bamfiles/$file\_ctrl.selected.bam -totalReads $norm 

done


#- (peak calling)
mkdir Peaks
blackList="/home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed"

for file in MLL1 MEN1 PolII-RPB1
do

	findPeaks TagDir/$file\_expr -i TagDir/$file\_ctrl -style histone -C 0 -o Peaks/$file\_broad 
	findPeaks TagDir/$file\_expr -i TagDir/$file\_ctrl -style factor -C 0 -o Peaks/$file\_narrow
	pos2bed.pl -o Peaks/$file\_broad.bed Peaks/$file\_broad
	pos2bed.pl -o Peaks/$file\_narrow.bed Peaks/$file\_narrow
	intersectBed -v -a Peaks/$file\_broad.bed -b $blackList | grep ^chr | grep -v chrM > Peaks/$file\_broad.Clean.bed
	intersectBed -v -a Peaks/$file\_narrow.bed -b $blackList | grep ^chr | grep -v chrM > Peaks/$file\_narrow.Clean.bed
	cat Peaks/$file\_broad.Clean.bed Peaks/$file\_narrow.Clean.bed | sort -k1,1 -k2,2n | bedtools merge -i - > Peaks/$file\_all.Clean.bed

done

for file in JunD
do

	findPeaks TagDir/$file\_expr -i TagDir/$file\_ctrl -style factor -C 0 -o Peaks/$file\_narrow
	pos2bed.pl -o Peaks/$file\_broad.bed Peaks/$file\_broad
	intersectBed -v -a Peaks/$file\_narrow.bed -b $blackList | grep ^chr | grep -v chrM > Peaks/$file\_narrow.Clean.bed
done 

#-- Overlapping peaks 
#__________________________________


Rscript ann.R -m 100 -d Peaks/h3k4me3.ENCFF862LUQ.replicatedNarrowPeaks.bed,./Peaks/JunD_narrow.Clean.bed,./Peaks/MLL1_all.Clean.bed -c 30 -n h3k4me3,JunD,MLL1 -o Overlap.bed -i ./Peaks/MEN1_all.Clean.bed

cat Overlap.bed.totalMatrix.txt | awk '{if($6==0 && $7==0 && $8==0) print $0}' | cut -f1-3 > ./Peaks/MEN1.unknown.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==1 && $7==1 && $8==1) print $0}' | cut -f1-3 > ./Peaks/MEN1.h3k4me3JunDMLL1.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==0 && $7==1 && $8==0) print $0}' | cut -f1-3 > ./Peaks/MEN1.onlyJunD.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==0 && $7==0 && $8==1) print $0}' | cut -f1-3 > ./Peaks/MEN1.onlyMLL1.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==1 && $7==0 && $8==0) print $0}' | cut -f1-3 > ./Peaks/MEN1.onlyh3k4me3.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==0 && $7==1 && $8==1) print $0}' | cut -f1-3 > ./Peaks/MEN1.JunDMLL1.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==1 && $7==0 && $8==1) print $0}' | cut -f1-3 > ./Peaks/MEN1.h3k4me3MLL1.bed 
cat Overlap.bed.totalMatrix.txt | awk '{if($6==1 && $7==1 && $8==0) print $0}' | cut -f1-3 > ./Peaks/MEN1.h3k4me3JunD.bed 


#-- Motif identification
#__________________________________
mkdir Motifs

echo "category,total,AP1,ATF,Unknown" > Motifs/Counts.txt

for file in unknown  h3k4me3JunDMLL1 onlyJunD onlyMLL1 onlyh3k4me3 JunDMLL1 h3k4me3MLL1 h3k4me3JunD
do

	echo "======================================================"
	total=$(wc -l Peaks/MEN1.$file.bed | cut -f1 -d " ")
	findMotifsGenome.pl Peaks/MEN1.$file.bed hg38 ./Motifs -size 400 -nomotif -p 30 -find Motifs/ap1.txt > ./Motifs/MEN1.$file.AP1.txt ; 
	ap1=$(cat ./Motifs/MEN1.$file.AP1.txt  | sed 1d | cut -f1 | sort | uniq | wc -l )
	findMotifsGenome.pl Peaks/MEN1.$file.bed hg38 ./Motifs -size 400 -nomotif -p 30 -find Motifs/atf7.txt > ./Motifs/MEN1.$file.ATF.txt ; 
	atf=$(cat ./Motifs/MEN1.$file.ATF.txt  | sed 1d | cut -f1 | sort | uniq | wc -l )
	unknown=$(echo "$total-($ap1+$atf)" | bc -l )
	echo "$file,$total,$ap1,$atf,$unknown" >> Motifs/Counts.txt

done 

#- pie chart
R

d= read.table("Motifs/Counts.txt")
jpeg("pie.jpeg",unit="in",res=300,height=2*3,width=2*2)
opar<-par(mfrow=c(3,2))

pie(as.numeric(d[2,-c(1,2)]),col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","Unk"),main=expression("Cluster 1"))
pie(as.numeric(d[6,-c(1,2)]),col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","Unk"),main=expression("Cluster 2")) 
pie(as.numeric(d[8,-c(1,2)]),col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","Unk"),main=expression("Cluster 3")) 
pie(as.numeric(d[3,-c(1,2)]),col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","Unk"),main=expression("Cluster 4")) 
pie(as.numeric(d[4,-c(1,2)]) + as.numeric(d[5,-c(1,2)]) + as.numeric(d[7,-c(1,2)]), col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","NA"),main=expression("Cluster 5")) 
pie(as.numeric(d[1,-c(1,2)]),col=c("#4f81bd","#19477f","#8dc1ff"),labels=c("AP1","ATF","Unk"),main=expression("Cluster 6"))
dev.off()

#-- HeatMap
#__________________________________

mkdir Heatmaps
mkdir bamcompare

for file in MEN1 MLL1 JunD PolII-RPB1
do
	h_ex=$(fgrep read1 Bamfiles/$file\_expr.flagstat.txt | cut -d " " -f1)
	h_ct=$(fgrep read1 Bamfiles/$file\_ctrl.flagstat.txt | cut -d " " -f1)
	s_ex=$(fgrep read1 SpikeIn/$file\_expr.flagstat.txt  | cut -d " " -f1)
	s_ct=$(fgrep read1 SpikeIn/$file\_ctrl.flagstat.txt  | cut -d " " -f1)
	norm=$(echo "1/(($s_ct/$h_ct)/($s_ex/$h_ex))" | bc -l)
	bamCompare -b1 Bamfiles/$file\_expr.bam -b2 Bamfiles/$file\_ctrl.bam -o bamcompare/$file.bw -p 30 -bl /home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed --scaleFactors 1:$norm --effectiveGenomeSize 2913022398

done

# In the case of histone marks (ENCODE) dataset, all default parameters are used to generate bamcompare files (which means without SpikeIn normalization or without  --scaleFactors)
#-- heatmap

bedfiles="Peaks/MEN1.h3k4me3JunDMLL1.bed Peaks/MEN1.JunDMLL1.bed Peaks/MEN1.h3k4me3JunD.bed Peaks/MEN1.onlyJunD.bed Peaks/MEN1.onlyMLL1.bed Peaks/MEN1.onlyh3k4me3.bed Peaks/MEN1.h3k4me3MLL1.bed Peaks/MEN1.unknown.bed"
samples="./bamcompare/MEN1.bw ./bamcompare/MLL1.bw ./bamcompare/JunD.bw ./bamcompare/ENCFF063XTI.rep1.filtered.bamcompare.bw ./bamcompare/ENCFF241VRU.rep1.filtered.bamcompare.bw ./bamcompare/ENCFF617YCQ.rep1.filtered.bamcompare.bw ./bamcompare/ENCFF113QJM.rep1.filtered.bamcompare.bw ./bamcompare/ENCFF725TAB.rep1.filtered.bamcompare.bw ./bamcompare/PolII-RPB1.bw"
output="Heatmaps/out"
name="MEN1 MLL1 JunD H3k4me3 H3k4me2 H3k4me1 H3k27ac H3k9ac PolII(RPB1)"

computeMatrix reference-point -R $bedfiles -b 5000 -a 5000 --missingDataAsZero -bl /home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed --smartLabels -p 35  -S $samples -o $output.gz --samplesLabel $name
plotHeatmap -m $output.gz -o $output.jpeg --colorMap RdBu --dpi 300 --refPointLabel "" --zMin -2 --zMax 5 --whatToShow "heatmap and colorbar"

#-- heatmap (atacSeq) (since ATAC seq do not have any control bamcoverage file was used which generated by default parameters)

bedfiles="Peaks/MEN1.h3k4me3JunDMLL1.bed Peaks/MEN1.JunDMLL1.bed Peaks/MEN1.h3k4me3JunD.bed Peaks/MEN1.onlyJunD.bed Peaks/MEN1.onlyMLL1.bed Peaks/MEN1.onlyh3k4me3.bed Peaks/MEN1.h3k4me3MLL1.bed Peaks/MEN1.unknown.bed"
samples="./bamcompare/SRR8171284.deduplicate.bamcoverage.bw"
output="Heatmaps/outATAC"
name="ATAC"

computeMatrix reference-point -R $bedfiles -b 5000 -a 5000 --missingDataAsZero -bl /home/sheikh/Databases/EncodeBlackListedRegions/hg38-blacklist.v2.bed --smartLabels -p 35  -S $samples -o $output.gz --samplesLabel $name
plotHeatmap -m $output.gz -o $output.jpeg --colorMap RdBu --dpi 300 --refPointLabel "" --zMin -2 --zMax 5 --whatToShow "heatmap and colorbar"

#-- Coverage around peaks of menin for MEN1 and menin mutations
#_______________________________________

#- tagDirectory of the selected reads ()
Samples	Spike	TotalReads	SelectedReads	SpikeIn/TotalReads	Ratio	HomerReads(for --totalReads)
MEN1	439997	7382489	7382489	0.059600089	0.768656573	9604404.959
MEN1-E255K	873721	18283880	7381148	0.047786411	0.958682589	7699261.552
MEN1-E359K	1243376	24272765	7381272	0.051225149	0.894326342	8253443.571
MEN1-E408Q	1100307	24017755	7381272	0.045812233	0.999994903	7381309.62
MEN1-R52G	1820594	31667197	7379316	0.057491479	0.796848517	9260625.881
#-
makeTagDirectory TagDir/MEN1-E255K_Selected Bamfiles/MEN1-E255K_expr.selected.bam  -totalReads 7699261.552 &
makeTagDirectory TagDir/MEN1-E359K_Selected Bamfiles/MEN1-E359K_expr.selected.bam  -totalReads 8253443.571 &
makeTagDirectory TagDir/MEN1-E408Q_Selected Bamfiles/MEN1-E408Q_expr.selected.bam  -totalReads 7381309.62 &
makeTagDirectory TagDir/MEN1-R52G_Selected Bamfiles/MEN1-R52G_expr.selected.bam  -totalReads  9260625.881 &
makeTagDirectory TagDir/MEN1_Selected Bamfiles/MEN1_expr.selected.bam  -totalReads 9604404.959 &

#-- differential peaks 
#__________________________________

mkdir DifferentialPeaks
cat Peaks/MEN1.onlyMLL1.bed Peaks/MEN1.onlyh3k4me3.bed Peaks/MEN1.h3k4me3MLL1.bed > temp.bed

for file in MEN1-E408Q MEN1-R52G MEN1-E255K MEN1-E359K
do
	getDifferentialPeaks Peaks/MEN1.onlyJunD.bed TagDir/MEN1_Selected TagDir/$file\_Selected -F 0 -P 1 -size 400 > DifferentialPeaks/$file-onlyJunD.txt
	getDifferentialPeaks temp.bed TagDir/MEN1_Selected TagDir/$file\_Selected -F 0 -P 1 -size 400  > DifferentialPeaks/$file-onlyMLL1Histone.txt
done 

R

files=c("MEN1-R52G","MEN1-E408Q","MEN1-E255K","MEN1-E359K")

for (i in c(1:length(files))){
	file=files[i]
	d=read.table(paste("DifferentialPeaks/",file,"-onlyJunD.txt",sep=""))
	d$color = ifelse(d$V10> 4 & d$V11 < 0.0001, "grey72", ifelse (d$V10 < -4 & d$V11 < 0.0001, "red4", "grey30") )
	m = d[d$color=="grey30",]
	n = d[d$color=="grey72",]
	write.table(m[,c(2,3,4)],paste("DifferentialPeaks/",file,"-junD_noloss.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
	write.table(n[,c(2,3,4)],paste("DifferentialPeaks/",file,"-junD_loss.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
	#--
	d=read.table(paste("DifferentialPeaks/",file,"-onlyMLL1Histone.txt",sep=""))
	d$color = ifelse(d$V10> 4 & d$V11 < 0.0001, "grey72", ifelse (d$V10 < -4 & d$V11 < 0.0001, "red4", "grey30") )
	m = d[d$color!="grey72",]
	n = d[d$color=="grey72",]
	write.table(m[,c(2,3,4)],paste("DifferentialPeaks/",file,"-onlyMLL1Histone_noloss.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
	write.table(n[,c(2,3,4)],paste("DifferentialPeaks/",file,"-onlyMLL1Histone_loss.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}

q("no")

#- barplot 
#_______________________________
Number of the peaks
#- 3201
   677 MEN1-E255K-onlyMLL1Histone_noloss.txt
   754 MEN1-E359K-onlyMLL1Histone_noloss.txt
  3146 MEN1-E408Q-onlyMLL1Histone_noloss.txt
   776 MEN1-R52G-onlyMLL1Histone_noloss.txt

#- 1187
   546 MEN1-E255K-junD_noloss.txt
   479 MEN1-E359K-junD_noloss.txt
  1140 MEN1-E408Q-junD_noloss.txt
   321 MEN1-R52G-junD_noloss.txt

R

histone=100*(c(3146,776,677,754)/3201)
jun=100*(c(1140,321,546,479)/1187)
d=t(as.matrix(data.frame(histone,jun)))
colnames(d) <- c("E408Q","R52G","E255K","E359K")
jpeg("barplotData-2Feb2022.jpeg",unit="in",res=300,height=4,width=4)
plot(c(0.5,12),c(0,110),type="n",axes=F,xlab="", ylab= "Percentage of peaks retained")
grid()
box()
x <- barplot(d, beside = T, las=2, ylim = c(0, 110), col = c("#c00000","#4f81bd"), add =T, border = "black", cex.axis = 1, cex = 1 )
legend("topright",col="black",pt.bg=c("#4f81bd","#c00000"), legend=c("Cluster 4","Cluster 5"),bty='n', pch = 22,cex = 1)

dev.off()
q("no")


