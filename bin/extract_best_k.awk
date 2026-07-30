#!/usr/bin/awk -f

/<h2>Predicted best k: [0-9]+<\/h2>/ {
    gsub(/.*Predicted best k: /, "")
    gsub(/<\/h2>.*/, "")
    print
    exit
}
