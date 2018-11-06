NPIC=1 # Change this number 1--6 (too lazy to implement loop in bash)
convert orig/greyFace${NPIC}_V3-01.jpg -sigmoidal-contrast 20,60% -resize 256x256 ${NPIC}.jpg
