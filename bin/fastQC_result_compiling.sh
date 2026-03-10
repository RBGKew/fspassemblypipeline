for i in `ls -d */`;do grep 'Filename' $i/*fastqc_data.txt ;done |awk '{split($2,a,".");print a[1]}'> name.txt
for i in `ls -d */`;do grep '%GC' $i/*fastqc_data.txt ;done |awk '{print $2}'> GC.txt
for i in `ls -d */`;do grep 'Total Sequences' $i/*fastqc_data.txt ;done |awk '{print $3}' > total_reads_no.txt
for i in `ls -d */`;do grep 'Total Deduplicated Percentage' $i/*fastqc_data.txt ;done  |awk '{print 100-$4}'> duplication_level.txt
for i in `ls -d */`;do grep 'Per sequence GC content' $i/*gz_summary.txt ;done|awk '{print $1}' > GC_pass.txt
echo "name	GC	GC_pass	total_reads_no	Duplicated Percentage" > fastQC_result.txt
paste name.txt GC.txt GC_pass.txt total_reads_no.txt duplication_level.txt >> fastQC_result.txt
rm name.txt GC.txt GC_pass.txt total_reads_no.txt duplication_level.txt
