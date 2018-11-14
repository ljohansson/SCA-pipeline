#!/bin/sh
for f in *.txt; do
    d="$(head -1 "$f" | awk '{print $1}').txt"
    if [ ! -f "$d" ]; then
        mv "$f" "$d"
    else
        echo "File '$d' already exists! Skiped '$f'"
    fi
    tail -n +2 $f > "$f.tmp" && mv "$f.tmp" "$f"
done

