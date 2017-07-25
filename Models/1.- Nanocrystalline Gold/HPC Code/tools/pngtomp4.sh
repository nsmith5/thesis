#!/bin/bash

if [ "$2" = "" ]; then
	output="out.mp4"
else
	output="$2"
fi

input="'$1/*.png'"

echo 'ffmpeg -r 10 -pattern_type glob -i $input -c:v libx264 -pix_fmt yuv420p $output'

ffmpeg -r 10 -pattern_type glob -i $input -c:v libx264 -pix_fmt yuv420p $output
