#!/bin/bash

mkdir -p ../build/doc/
cp -r img/ ../build/doc/

pandoc -s -S -f markdown -t html5 --mathjax --toc \
    --biblatex --bibliography=./refs.bib PclePush.md > test.html
