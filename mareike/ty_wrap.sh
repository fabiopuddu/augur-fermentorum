# /bin/sh

# written 22.09.2016 


for fol in Del*_*

do date >> log

echo $fol >> log

cd $fol 

/mnt/home3/jackson/fp305/sw/bin/PF/shell/ty-realign.sh

sleep 8

cd ..

done

