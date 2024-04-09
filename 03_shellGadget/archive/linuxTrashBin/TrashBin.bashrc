alias rm='trash_func'
# Alias: rm
# Description: Move files to trash
# eg: rm file1.txt file2.jpg

alias rr='retrieve_func'
# Alias: rr
# Description: Restore files or directories from trash
# eg: rr file.txt

alias rc='cleartrash_func'
# Alias: rc
# Description: Clear all files and directories from trash

alias rs='searchtrash_func'
# Alias: rs
# Description: Search for files in trash, can optionally display file details
# eg: rs file.txt; rs file.txt --info; rs --all --info

alias rd='searchtrash_bydate_func'
# Alias: rd
# Description: Search for files in trash by date
# eg: rd 2022-01-01

function trash_func() {
    mkdir -p ~/.trash  # Create the trash directory

    files=("$@")
    for i in "${files[@]}"; do
        if [[ -e $i ]]; then
            if [[ -e ~/.trash/$i ]]; then
                /bin/rm -r ~/.trash/$i
            fi
            cp -Rp $i ~/.trash/$(date +%m-%d-%H-%M-%S)-$(basename $i)
            mv $i ~/.trash/ || echo "Failed to move $i to trash!"
        else
            echo "warning: file >>> $i <<< not found!"
        fi
    done
}

function retrieve_func() {
    if [[ -z $1 ]]; then
        echo "Please provide the filename or directory name to restore"
        return 0
    fi

    if [[ $1 == "all" ]]; then
        mv -i ~/.trash/* ./ || echo "Failed to retrieve files from trash!"
    else
        files=("$@")
        for i in "${files[@]}"; do
                # Restore file/folder
                mv -i ~/.trash/$i ./ || echo "Failed to retrieve $i from trash!"
        done
    fi
}

function cleartrash_func() {
    echo -n "Y/y: "
    read answer
    if [[ $answer == [Yy] ]]; then
        /bin/rm -rf ~/.trash/*
    else
        echo "Operation failed"
    fi
}


function searchtrash_func() {
    if [[ -z $1 ]]; then
        echo "Please provide the file name to search for"
        return 0
    fi

    if [[ $1 == "--all" ]]; then
        found_files=$(find ~/.trash -type f)
    else
        search_term=$1
        found_files=$(find ~/.trash -iname "*$search_term*")
    fi

    if [[ -z $found_files ]]; then
        echo "No files found containing '$1'"
    else
        echo "Found files:"
        echo "$found_files"

        if [[ $2 == "--info" ]]; then
            while IFS= read -r file_path; do
                echo
                echo "File details:"
                ls -l "$file_path"
                echo
                file "$file_path"
            done <<< "$found_files"
        fi
    fi
}


function searchtrash_bydate_func() {
    if [[ -z $1 ]]; then
        echo "Please provide the date in the format 'YYYY-MM-DD' to search for"
        return 0
    fi

    date_to_search=$1
    if ! [[ $date_to_search =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
        echo "Invalid date format. Please provide the date in the format 'YYYY-MM-DD'"
        return 0
    fi

    found_files=$(find ~/.trash -type f -newermt $date_to_search -exec sh -c 'if [ "$(date -r "$1" +%Y-%m-%d)" = $2 ]; then echo "$1"; fi' sh {} $date_to_search \;)
    if [[ -z $found_files ]]; then
        echo "No files found created on '$date_to_search'"
    else
        echo "Found files created on '$date_to_search':"
        echo "$found_files"
    fi
}