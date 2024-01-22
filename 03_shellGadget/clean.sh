#!/bin/bash
#Desc: 清理文件夹和文件
#Params: -d: 只清理文件夹 -f: 只清理文件，-s: 后面的的第一个参数接的是后缀名，只清理指定后缀文件，-i: 后面的第一个参数接的是文件中包含的字符，只清除包含指定字符的文件或者文件夹。剩余的参数是排除之外的文件或者文件夹
#Author: gchen

# 默认清理文件和文件夹
clean_files=false
clean_dirs=false

# 是否清理指定后缀
clean_dest_suffix=false

# 模糊清理
fuzzy_clean=false

# 解析命令行参数
args=`getopt -o afds:i: -- "$@"`
eval set -- "$args"
while true; do
  case ${1} in
    -a)
      clean_files=true
      clean_dirs=true
      shift 1
      ;;
    -f)
      clean_files=true
      shift 1
      ;;
    -d)
      clean_dirs=true
      shift 1
      ;;
    -s)
      clean_dest_suffix=true
      suffix=$2
      shift 2
      ;;
    -i)
      fuzzy_clean=true
      include_str=$2
      shift 2
      ;;
    --)
      shift
      break
      ;;
     *)
      echo "无效的选项 $1" >&2
      exit 1
      ;;
  esac
done

# 随机数
ran=${RANDOM}${RANDOM}
mkdir ~/.trash/${ran}

# 清理文件夹
if [ "$clean_dirs" = true ]; then
  echo "开始清理文件夹..."
  for dir in *; do
    if [ -d "$dir" ]; then

      flag=0
      for file_ in $@; do
          if [ "$dir" = "$file_" ];then
              flag=1
              break
          fi 
      done

      if [ "$flag" -eq 0 ];then
          mv "$dir" ~/.trash/${ran}
          echo "已删除文件夹: $dir"
      fi

    fi
  done
fi

# 清理文件
if [ "$clean_files" = true ]; then
  echo "开始清理文件..."

  for file in *; do
    if [ -f "$file" ] && [[ "./$file" != "$0" ]];then

      flag=0
      for file_ in $@; do
          if [ "$file" = "$file_" ];then
              flag=1
              break
          fi 
      done

      if [ "$flag" -eq 0 ];then
          mv "$file" ~/.trash/${ran}
          echo "已删除文件: $file"
      fi

    fi
  done
fi


# 清理指定后缀的文件
if [ "$clean_dest_suffix" = true ]; then
  echo "开始清理文件..."

  for file in *;do
      if [ -f ${file} ]  && [[ ${file} == *\.${suffix} ]];then

          flag=0
          for file_ in $@; do
              if [ "$file" = "$file_" ];then
                  flag=1
                  break
              fi 
          done

          if [ "$flag" -eq 0 ];then
              mv "$file" ~/.trash/${ran}
              echo "已删除文件: $file"
          fi

      fi
  done

fi


# 模糊清理 
if [ "$fuzzy_clean" = true ]; then
  echo "开始清理文件..."

  for file in *;do
      if [[ ${file} == *${include_str}* ]];then

          flag=0
          for file_ in $@; do
              if [ "$file" = "$file_" ];then
                  flag=1
                  break
              fi 
          done

          if [ "$flag" -eq 0 ];then
              mv "$file" ~/.trash/${ran}
              echo "已删除文件: $file"
          fi

      fi
  done

fi
   

echo "清理完成！"
