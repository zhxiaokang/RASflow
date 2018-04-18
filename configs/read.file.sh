while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Text read from file: $line"
    done < "$1"
