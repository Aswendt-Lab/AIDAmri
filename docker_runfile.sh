#!/bin/bash
SRCPATH=${@: -1}
dockerrun()
{
if [ "$bflag" == 'b' ]; then
	docker build -t aidamri:latest -f Dockerfile .
fi
docker run \
-dit --rm \
--name aidamri \
--mount type=bind,source=$SRCPATH,target=/aida/mountdata \
aidamri:latest
if [ "$aflag" == 'a' ]; then
	docker attach aidamri
fi
}


while getopts ":ba" opt; do
  case $opt in
    b)
		bflag='b'
    	;;
	a)
		aflag='a'
		;;
    \?)
    	echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
  esac
done


dockerrun