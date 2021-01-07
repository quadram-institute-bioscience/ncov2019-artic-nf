#!/usr/bin/env bash
rm .nextflow.* 2> /dev/null
rm *.txt.{1..9} 2> /dev/null
rm *.html.{1..9} 2> /dev/null

echo "Cleaned!"