alias rm=trash
alias rr=retrieve
alias rl='ls ~/.trash'
alias rc=cleartrash
function trash()
{
    if [ ! -d ~/.trash ];then
        mkdir ~/.trash
    fi

    for i in "$@"
    do
        if test -e ${i};then
            if test -e ~/.trash/${i};then
                /bin/rm -r ~/.trash/${i}
            fi
            cp -r $i ~/.trash/`date +%m-%d-%H-%M-%S`-`basename $i`
            mv $i ~/.trash/
        else
            echo warning: file \>\>\> ${i} \<\<\< not find!
        fi
    done
}
function retrieve()
{
    if [ -z "$1" ];then
        echo "请输入需要还原的文件名"
        return 0
    fi

    for i in "$@"
    do
        mv -i ~/.trash/$i ./
    done
}

function cleartrash()
{
    read -p "Y/N: " answer
    if [ "${answer}" = 'Y' -o "${answer}" = 'y' ];then
        /bin/rm -rf ~/.trash/*
    else
        echo "operate withdraw"
    fi
}