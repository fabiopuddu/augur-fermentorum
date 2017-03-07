# /bin/sh

# written 22.09.2016 
for f in Del*_*
	do
		cd "$f"
		echo "Submitting $f"
		ty-realign.sh
		cd ..
		sleep 1
	done
