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
        ref1="$(random_DNA 100)"
        printf "%s\n" ${ref1} >&3
        marker="${ref1:10:20}"
        printf ">%d\n%s\n" "${g}" "${marker}" >&4
        ref2="$(random_DNA 100)"
        printf "%s\n" ${ref2} >&5
        ref1=${ref1:0:60}
        ref2=${ref2:40:100}
        printf "0\t%s\t50\t10\t%s\t60\n" "${ref1}" "${ref2}" >&6
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

generate_seq_groupby_marker 3> example.ref1 4> example.fa 5> example.ref2 6> example.ref | make_fastq > example.fq
