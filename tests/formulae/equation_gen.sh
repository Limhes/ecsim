#!/bin/sh

if [ ! $# == 2 ]; then
  echo "Usage: $0 latex_equation output_filename"
  exit
fi

TEX_EQUATION="$1"

pdflatex "\def\formula{${TEX_EQUATION}}\input{formula.tex}" # > /dev/null
convert -density 150 formula.pdf -quality 100 -strip "$2" # > /dev/null
