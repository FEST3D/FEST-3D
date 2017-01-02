#!bin/dir

#src=src/
#
#subdir_list=$(ls -d ../$src/*/)
#echo $subdir_list
#
#check_dir=$subdir_list
#function  create_dir{
#  for dir in ${check_dir[@]};
#    if [ ! -d ../$src$dir ]; then
#      mkdir ../$src$dir
#      new_subdir=$(ls -d ../$src$dir*/)
#
#    fi
#
#  done
#}

function recurse() {
  for i in "$1"/*;do
    if [[ "$1/*" = "$i" ]];then
      continue
    fi
    if [ -d "$i" ];then
      dir=${$1:6:}
      echo "dir: $dir"
      recurse "$i"
    elif [ -f "$i" ]; then
      #echo "file: $i"
      continue
    fi
  done
}
recurse ../src
