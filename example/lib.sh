random_DNA() {
    local length=$1
    local chars="ACGT"
    local str=""
    for ((i = 0; i < ${length}; ++i))
    do
        str+=${chars:RANDOM%${#chars}:1}
    done
    printf "$str"
}

make_fastq() {
    awk '
        {
            printf("@seq%d\n%s\n+\n", NR - 1, $0)
            for (i=0; i<length($0); ++i) {
                printf("~")
            }
            printf("\n")
        }
    '
}

generate_seq_groupby_marker() {
    for ((g = 0; g < 10; ++g))
    do
        marker="$(random_DNA 20)"
        printf ">%d\n%s\n" "${g}" "${marker}" >&3
        prefix="$(random_DNA 10)"
        suffix="$(random_DNA 70)"
        ref="${prefix}${marker}${suffix}"
        ref1=${ref:0:60}
        ref2=${ref:40:100}
        printf "0\t%s\t50\t10\t%s\t60\n" "${ref1}" "${ref2}" >&4
        for ((s = 0; s < 10; ++s))
        do
            up="$(( "${RANDOM}" % 21 - 10 ))"
            up_seq="${ref1:0:50+${up}}"
            down="$(( "${RANDOM}" % 21 - 10 ))"
            down_seq="${ref2:${#ref2}-50+${down}}"
            random_len="$(("${RANDOM}" % 5))"
            random_seq="$(random_DNA ${random_len})"
            printf "%s\n" "${up_seq}${random_seq}${down_seq}"
        done
    done
}

# 切换运行路径到脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

generate_seq_groupby_marker 3> example.fa 4> example.ref | make_fastq > example.fq
