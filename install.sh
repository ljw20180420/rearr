#!/bin/bash
# simply append three paths to ~/.bashrc

# echo "your conda environments are:"
# conda env list | sed '/^#/d' |sed '/^$/d' | cut -d" " -f1 | nl -w1
# env_list=$(conda env list | sed '/^#/d' |sed '/^$/d' | cut -d" " -f1)
# default_env=$(echo $env_list | cut -d" " -f1)
# read -p "which environment to use (you can also create a new one)? (default: $default_env):" chosen_env
# if [[ $chosen_env == "" ]]
# then
#     echo "you do not choose, use the default env $default_env"
#     chosen_env=$default_env
# else
#     echo "you choose $chosen_env"
#     if [[ ! $env_list =~ $chosen_env ]]
#     then
#         read -p "$chosen_env do not exist, do you want to create it? [y/n] (default: n):" create_env
#         case $create_env in
#             y|Y|yes|Yes|YES|yeah|Yeah|YEAH)
#                 conda create -n $chosen_env more-itertools biopython numpy;;
#         esac
#         exit
#     fi  
# fi
# conda install -n $chosen_env more-itertools biopython numpy

cd $(dirname $0)
if [[ $(grep -c $(pwd) ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi
if [[ $(grep -c $(pwd)/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
fi
cmake -DCMAKE_BUILD_TYPE=Release Rearrangement/build
make -C Rearrangement/build
if [[ $(grep -c $(pwd)/Rearrangement/build ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/Rearrangement/build\":"'$PATH' >> ~/.bashrc
fi