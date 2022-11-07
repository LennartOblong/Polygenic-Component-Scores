#Takes in the melodic_oIC.nii.gz file,
#Which is a 4D volume of all components (raw ICs, NOT z-transformed)
#And outputs text files with SNP IC values in rows
#Each text file is a single IC


input_melodic_oIC=$1
output_fname=$2
variant_base_file=variants-17103079-clump-masked.txt

n_ic=$(fslinfo $input_melodic_oIC |awk '$1=="dim4"{print $2}')

fslsplit $input_melodic_oIC split_tmp -t

for ((i=1;i<=n_ic;i++));do
fsl2ascii split_tmp$(printf "%04d" $((i-1)) ) split_tmp$(printf "%04d" $((i-1)) )"_ascii"

cat split_tmp$(printf "%04d" $((i-1)) )"_ascii00000"|tr "\n" " "|sed 's/   / /g'|sed 's/  / /g'|sed 's/ /\n/g'> split_tmp$(printf "%04d" $((i-1)) )"_ascii"

#only keep the in-mask SNPs in the output text
paste -d" " $variant_base_file split_tmp$(printf "%04d" $((i-1)) )"_ascii"| awk '$6==1{print $7}' > ${output_fname}_ic_$i
rm split_tmp$(printf "%04d" $((i-1)) )"_ascii" split_tmp$(printf "%04d" $((i-1)) )"_ascii00000" split_tmp$(printf "%04d" $((i-1)) )".nii.gz" 
done

